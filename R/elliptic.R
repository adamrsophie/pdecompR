#' Linear elliptic problems
#'
#' `elliptic_problem()` defines the scalar second-order problem
#'
#' `-div(k grad u) + b . grad u + c u = f` in the domain,
#'
#' with Dirichlet, Neumann and Robin boundary conditions. `poisson_problem()`
#' is the special case `b = 0`, `c = 0`.
#'
#' The weak form is: find `u` with the Dirichlet data such that for all test
#' functions `v` vanishing on the Dirichlet boundary
#'
#' `int k grad u . grad v + (b . grad u) v + c u v + sum_Robin int alpha u v`
#' `= int f v + sum_Neumann int g v + sum_Robin int r v`.
#'
#' Coefficients may be numbers or vectorized `function(x, y)`; `diffusion`
#' may also be a constant symmetric positive definite 2 by 2 matrix, and
#' `advection` a length-2 vector or a function returning a two-column matrix.
#' Uniform ellipticity (`k > 0`) is checked at the quadrature points.
#'
#' With `stabilization = "supg"` (P1 only) the streamline upwind
#' Petrov-Galerkin terms
#' `sum_T tau_T int_T (b . grad u + c u - f) (b . grad v)` are added, with
#' `tau_T = h_T / (2 |b|) (coth(Pe) - 1 / Pe)`, `Pe = |b| h_T / (2 k)`,
#' evaluated at the element centroid with `h_T` the longest edge. The method
#' is consistent: the exact solution satisfies the stabilized equations.
#'
#' Problems without Dirichlet or Robin conditions and without reaction are
#' pure Neumann problems: the solution is unique only up to a constant. The
#' solver then checks the compatibility condition `int f + int g = 0` and
#' returns the solution with zero mean, labeled as such. The discrete check is
#' `|sum_i b_i| <= compatibility_tol * sum_i |b_i|` (an argument of
#' [solve_pde()], default `1e-6`), which tolerates quadrature error in
#' compatible data; a Lagrange multiplier for the zero-mean constraint absorbs
#' the remaining defect and is reported.
#'
#' @param mesh A [pde_mesh()].
#' @param diffusion Diffusion coefficient `k`.
#' @param reaction Reaction coefficient `c` (non-negative for coercivity).
#' @param advection Velocity `b`, or `NULL`.
#' @param source Source `f`.
#' @param bcs A `pde_bc` or list of them. Boundary edges without a condition
#'   receive the natural condition `k du/dn = 0`.
#' @param boundary Alias of `bcs` for [poisson_problem()].
#' @param degree Lagrange degree (1 or 2).
#' @param stabilization `"none"` or `"supg"`.
#' @param quad_degree Quadrature exactness degree (default `2 degree + 2`).
#' @return An object of class `pde_problem` (subclass `elliptic_problem`).
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' prob <- poisson_problem(mesh, source = function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y),
#'                         boundary = dirichlet_bc(0))
#' sol <- solve_pde(prob)
#' sol
#' @export
elliptic_problem <- function(mesh, diffusion = 1, reaction = 0, advection = NULL,
                             source = 0, bcs = list(), degree = 1,
                             stabilization = c("none", "supg"), quad_degree = NULL) {
  stopifnot(inherits(mesh, "pde_mesh"))
  stabilization <- match.arg(stabilization)
  bcs <- normalize_bcs(bcs)
  check_bc_markers(bcs, mesh)
  degree <- check_count(degree, "degree")
  if (stabilization == "supg") {
    if (is.null(advection)) stopf("SUPG stabilization needs an advection field.")
    if (degree != 1L) stopf("SUPG stabilization is implemented for P1 elements only.")
  }
  check_coefficient(diffusion, "diffusion", allow_matrix = TRUE)
  check_coefficient(reaction, "reaction")
  check_coefficient(source, "source")
  structure(list(mesh = mesh, diffusion = diffusion, reaction = reaction,
                 advection = advection, source = source, bcs = bcs, degree = degree,
                 stabilization = stabilization, quad_degree = quad_degree),
            class = c("elliptic_problem", "pde_problem"))
}

#' @rdname elliptic_problem
#' @export
poisson_problem <- function(mesh, diffusion = 1, source = 0, boundary = list(), degree = 1,
                            bcs = boundary) {
  elliptic_problem(mesh, diffusion = diffusion, source = source, bcs = bcs, degree = degree)
}

check_coefficient <- function(v, name, allow_matrix = FALSE) {
  if (is.null(v) || is.function(v) || inherits(v, "fe_function")) return(invisible(TRUE))
  if (allow_matrix && is.matrix(v) && is.numeric(v) && all(dim(v) == 2L)) return(invisible(TRUE))
  if (is.numeric(v) && length(v) == 1L && is.finite(v)) return(invisible(TRUE))
  stopf("`%s` must be a finite number or a function of (x, y).", name)
}

#' @export
print.pde_problem <- function(x, ...) {
  cat("<", class(x)[1], "> on a mesh with ", nrow(x$mesh$triangles),
      " triangles, P", x$degree, "\n", sep = "")
  if (length(x$bcs)) {
    for (bc in x$bcs) cat("  ", bc$type, " on ", paste(bc$marker, collapse = ", "), "\n", sep = "")
  } else {
    cat("  natural (zero flux) boundary conditions everywhere\n")
  }
  invisible(x)
}

#' Replace the mesh of a problem
#'
#' Returns a copy of `problem` posed on `mesh`; coefficients and boundary
#' conditions are kept. Used by convergence studies and adaptive loops.
#'
#' @param problem A `pde_problem`.
#' @param mesh The new [pde_mesh()].
#' @return The updated problem.
#' @export
update_mesh <- function(problem, mesh) {
  stopifnot(inherits(problem, "pde_problem"), inherits(mesh, "pde_mesh"))
  check_bc_markers(problem$bcs, mesh)
  problem$mesh <- mesh
  problem
}

# Element-wise SUPG parameter (P1).
supg_tau <- function(mesh, diffusion, advection, degree = 1L) {
  P <- mesh$nodes; T <- mesh$triangles
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
  q <- mesh_quality(mesh)
  h <- q$h_max
  if (is.numeric(advection)) {
    bn <- rep(sqrt(sum(advection^2)), length(cx))
  } else {
    bv <- as.matrix(advection(cx, cy)); bn <- sqrt(rowSums(bv^2))
  }
  k <- if (is.matrix(diffusion)) rep(max(eigen(diffusion, only.values = TRUE)$values), length(cx))
       else eval_spatial(diffusion, cx, cy, name = "diffusion")
  pe <- bn * h / (2 * degree * k)
  xi <- ifelse(pe > 1e-6, 1 / tanh(pe) - 1 / pe, pe / 3)
  ifelse(bn > 0, h / (2 * degree * pmax(bn, .Machine$double.xmin)) * xi, 0)
}

#' Discretize a PDE problem
#'
#' Builds the finite-element space and assembles the global system of a
#' linear stationary problem. The result separates the numerical
#' representation from the mathematical problem definition.
#'
#' @param problem A `pde_problem`.
#' @param ... Unused.
#' @return An object of class `pde_discretization` with the space, the full
#'   matrix `A` and load `b` (before Dirichlet elimination), the Dirichlet
#'   degrees of freedom and values, and assembly timings.
#' @export
discretize <- function(problem, ...) UseMethod("discretize")

#' @export
discretize.elliptic_problem <- function(problem, t = NULL, ...) {
  t0 <- proc.time()[["elapsed"]]
  V <- fe_space(problem$mesh, problem$degree)
  geo <- element_geometry(V, problem$quad_degree)
  D <- diffusion_qp(problem$diffusion, geo, t)
  kmin <- if (is.list(D)) NA else min(D)
  if (!is.list(D) && kmin <= 0) {
    stopf("The diffusion coefficient must be positive (minimum %s at quadrature points).",
          fmt_num(kmin))
  }
  cq <- coef_qp(problem$reaction, geo, t, "reaction")
  if (all(cq == 0)) cq <- NULL
  bq <- vec_qp(problem$advection, geo, t)
  fq <- coef_qp(problem$source, geo, t, "source")
  advt <- NULL; gload <- NULL
  if (problem$stabilization == "supg") {
    tau <- supg_tau(problem$mesh, problem$diffusion, problem$advection)
    tq <- matrix(tau, geo$m, geo$nq)
    D <- list(xx = D + tq * bq$x * bq$x, xy = tq * bq$x * bq$y,
              yx = tq * bq$y * bq$x, yy = D + tq * bq$y * bq$y)
    if (!is.null(cq)) advt <- list(x = tq * cq * bq$x, y = tq * cq * bq$y)
    gload <- list(x = tq * fq * bq$x, y = tq * fq * bq$y)
  }
  A <- global_matrix(V, local_bilinear(geo, mass = cq, adv = bq, advt = advt, diff = D))
  b <- global_vector(V, local_linear(geo, f = fq, g = gload))
  nat <- natural_bc_terms(V, problem$bcs, t, problem$quad_degree)
  if (!is.null(nat$A)) A <- A + nat$A
  b <- b + nat$b
  dd <- dirichlet_data(V, problem$bcs, t)
  has_robin <- any(vapply(problem$bcs, function(bc) bc$type == "robin", logical(1)))
  structure(list(space = V, A = A, b = b, dirichlet = dd,
                 free = setdiff(seq_len(V$ndof), dd$dofs),
                 pure_neumann = !length(dd$dofs) && !has_robin && is.null(cq),
                 assembly_time = proc.time()[["elapsed"]] - t0),
            class = "pde_discretization")
}

#' @export
print.pde_discretization <- function(x, ...) {
  cat("<pde_discretization> P", x$space$degree, ", ", x$space$ndof, " dofs (",
      length(x$free), " free), nnz(A) = ", length(x$A@x), "\n", sep = "")
  invisible(x)
}

#' Solve a PDE problem
#'
#' Generic solver entry point. Methods exist for elliptic, variational,
#' parabolic, reaction-diffusion and nonlinear problems.
#'
#' @param problem A `pde_problem`.
#' @param ... Passed to methods: `method`, `tol` and other arguments of
#'   [solve_linear_system()] for linear problems, time stepping or Newton
#'   controls for the others.
#' @return A `pde_solution` (stationary problems) or `pde_transient`
#'   (time-dependent problems).
#' @export
solve_pde <- function(problem, ...) UseMethod("solve_pde")

#' @export
solve_pde.elliptic_problem <- function(problem, method = "auto", ...) {
  disc <- discretize(problem)
  sol <- solve_discretization(disc, method = method, ...)
  sol$problem <- problem
  sol
}

# Dirichlet elimination and linear solve shared by all linear stationary
# problems.
solve_discretization <- function(disc, method = "auto", compatibility_tol = 1e-6, ...) {
  V <- disc$space
  u <- numeric(V$ndof)
  dd <- disc$dirichlet
  u[dd$dofs] <- dd$values
  free <- disc$free
  A <- disc$A
  t0 <- proc.time()[["elapsed"]]
  gauge <- NULL
  if (disc$pure_neumann) {
    one <- assemble_vector(V, 1)
    comp <- sum(disc$b)
    scale <- sum(abs(disc$b))
    if (abs(comp) > compatibility_tol * max(scale, 1e-300)) {
      stopf(paste0("Pure Neumann problem with incompatible data: ",
                   "int f + int g = %s (relative %s). The problem has no solution."),
            fmt_num(comp), fmt_num(comp / max(scale, 1e-300)))
    }
    K <- rbind(cbind(A, one), c(one, 0))
    s <- solve_linear_system(K, c(disc$b, 0), method = "lu", ...)
    u <- s$x[seq_len(V$ndof)]
    gauge <- list(type = "zero mean", multiplier = s$x[V$ndof + 1L],
                  compatibility = comp)
    info <- s$info
  } else {
    if (!length(free)) {
      info <- list(method = "none", status = "solved", residual = 0, relative_residual = 0,
                   iterations = 0L, n = 0L, nnz = 0L, symmetric = TRUE, tol = NA)
    } else {
      Aff <- A[free, free, drop = FALSE]
      rhs <- disc$b[free]
      if (length(dd$dofs)) rhs <- rhs - as.numeric(A[free, dd$dofs, drop = FALSE] %*% dd$values)
      s <- solve_linear_system(Aff, rhs, method = method, ...)
      u[free] <- s$x
      info <- s$info
    }
  }
  solve_time <- proc.time()[["elapsed"]] - t0
  diag <- structure(list(
    ndof = V$ndof, free = length(free), nnz = length(A@x),
    assembly_time = disc$assembly_time, solve_time = solve_time,
    solver = info, gauge = gauge,
    full_residual = if (length(free)) sqrt(sum((disc$b[free] - as.numeric(A[free, , drop = FALSE] %*% u))^2)) else 0
  ), class = "pde_diagnostics")
  structure(list(u = fe_function(V, u), space = V, diagnostics = diag, discretization = disc),
            class = "pde_solution")
}

#' @export
print.pde_solution <- function(x, ...) {
  d <- x$diagnostics
  cat("<pde_solution> P", x$space$degree, ", ", d$ndof, " dofs; solver ", d$solver$method,
      " (", d$solver$status, "), relative residual ", fmt_num(d$solver$relative_residual, 3),
      "\n", sep = "")
  if (!is.null(d$gauge)) cat("  pure Neumann problem: solution normalized to zero mean\n")
  cat("  u range [", fmt_num(min(x$u$values)), ", ", fmt_num(max(x$u$values)), "]\n", sep = "")
  invisible(x)
}

#' @export
summary.pde_solution <- function(object, ...) {
  d <- object$diagnostics
  structure(list(degree = object$space$degree, mesh = summary(object$space$mesh),
                 diagnostics = d, range = range(object$u$values)),
            class = "summary.pde_solution")
}

#' @export
print.summary.pde_solution <- function(x, ...) {
  d <- x$diagnostics
  cat("Finite-element solution (Lagrange P", x$degree, ")\n", sep = "")
  cat("  mesh: ", x$mesh$nodes, " nodes, ", x$mesh$elements, " triangles, h = ",
      fmt_num(x$mesh$h_max), "\n", sep = "")
  cat("  degrees of freedom: ", d$ndof, " (", d$free, " free)\n", sep = "")
  cat("  matrix nonzeros: ", d$nnz, "\n", sep = "")
  cat("  solver: ", d$solver$method, ", status ", d$solver$status, sep = "")
  if (!is.na(d$solver$iterations)) cat(", ", d$solver$iterations, " iterations", sep = "")
  cat("\n  linear residual: ", fmt_num(d$solver$residual, 3), " (relative ",
      fmt_num(d$solver$relative_residual, 3), ")\n", sep = "")
  if (!is.null(d$solver$note)) cat("  note: ", d$solver$note, "\n", sep = "")
  if (!is.null(d$gauge)) {
    cat("  pure Neumann: zero-mean gauge, compatibility defect ",
        fmt_num(d$gauge$compatibility, 3), "\n", sep = "")
  }
  cat("  assembly ", fmt_num(d$assembly_time, 3), " s, solve ", fmt_num(d$solve_time, 3),
      " s\n", sep = "")
  cat("  solution range [", fmt_num(x$range[1]), ", ", fmt_num(x$range[2]), "]\n", sep = "")
  cat("  (the linear residual measures the algebraic solve, not the discretization error)\n")
  invisible(x)
}
