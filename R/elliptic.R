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
#' evaluated at the element centroid with `h_T` the longest edge (for a
#' matrix `K`, `k` is the streamline diffusivity `b'Kb / |b|^2`). For P1
#' elements the diffusive part of the element residual vanishes, so the
#' method is consistent when the diffusion is constant; with variable
#' diffusion the term `grad k . grad u_h` is omitted, which keeps the
#' convergence order but is not strictly consistent.
#'
#' When no Dirichlet data are imposed and constants lie in the kernel of the
#' assembled operator (no reaction, no Robin term with positive `alpha`), the
#' problem is a pure Neumann problem: the solution is unique only up to a
#' constant. The solver then checks the compatibility condition (for a
#' symmetric operator `int f + int g = 0`) and
#' returns the solution with zero mean, labeled as such. The compatibility
#' defect is read from the Lagrange multiplier of the zero-mean constraint and
#' compared with `compatibility_tol * sum_i |b_i|` (an argument of
#' [solve_pde()], default `1e-6`). For a nonsymmetric operator (advection) the
#' left null vector is not constant; there a defect above `h^2` gives a
#' warning instead of an error, because data that are compatible for the PDE
#' leave a defect of the size of the discretization error; defects above
#' about one percent are still errors.
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
#' @references Brooks, A. N. and Hughes, T. J. R. (1982). Streamline
#'   upwind/Petrov-Galerkin formulations for convection dominated flows.
#'   *Computer Methods in Applied Mechanics and Engineering*, 32, 199-259.
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
  if (allow_matrix && inherits(v, "pde_tensor")) return(invisible(TRUE))
  if (is.numeric(v) && length(v) == 1L && is.finite(v)) return(invisible(TRUE))
  stopf("`%s` must be a finite number or a function of (x, y).", name)
}

#' @export
print.pde_problem <- function(x, ...) {
  cat("<", class(x)[1], "> on a mesh with ", nrow(x$mesh$triangles),
      " triangles, P", x$degree, "\n", sep = "")
  show <- function(bcs, prefix = "") {
    for (bc in bcs) cat("  ", prefix, bc$type, " on ", paste(bc$marker, collapse = ", "), "\n", sep = "")
  }
  if (inherits(x, "reaction_diffusion_problem")) {
    for (s in x$species) show(x$bcs[[s]], paste0(s, ": "))
    if (!length(unlist(x$bcs, recursive = FALSE))) cat("  natural (zero flux) boundary conditions everywhere\n")
  } else if (length(x$bcs)) {
    show(x$bcs)
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
#' @examples
#' prob <- poisson_problem(pde_mesh_rectangle(nx = 2), source = 1, boundary = dirichlet_bc(0))
#' finer <- update_mesh(prob, pde_mesh_rectangle(nx = 8))
#' finer
#' @export
update_mesh <- function(problem, mesh) {
  stopifnot(inherits(problem, "pde_problem"), inherits(mesh, "pde_mesh"))
  if (inherits(problem, "variational_problem")) {
    stopf(paste0("The forms of a variational problem are tied to their function space; ",
                 "rebuild them on the new mesh instead of calling update_mesh()."))
  }
  if (inherits(problem, "reaction_diffusion_problem")) {
    for (b in problem$bcs) check_bc_markers(b, mesh)
  } else {
    check_bc_markers(problem$bcs, mesh)
  }
  problem$mesh <- mesh
  problem
}

# Diffusion at quadrature points plus the SUPG streamline term tau b b'.
supg_diffusion <- function(D, tq, bq) {
  if (!is.list(D)) D <- list(xx = D, xy = 0, yx = 0, yy = D)
  list(xx = D$xx + tq * bq$x * bq$x, xy = D$xy + tq * bq$x * bq$y,
       yx = D$yx + tq * bq$y * bq$x, yy = D$yy + tq * bq$y * bq$y)
}

# Element-wise SUPG parameter (P1).
supg_tau <- function(mesh, diffusion, advection, degree = 1L) {
  P <- mesh$nodes; T <- mesh$triangles
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
  q <- mesh_quality(mesh)
  h <- q$h_max
  bv <- if (is.numeric(advection)) matrix(advection, length(cx), 2, byrow = TRUE)
        else as.matrix(advection(cx, cy))
  bn <- sqrt(rowSums(bv^2))
  K <- tensor_at(diffusion, cx, cy)
  k <- if (!is.null(K)) {
    # streamline diffusivity b'Kb / |b|^2
    sk <- (K$xx * bv[, 1]^2 + (K$xy + K$yx) * bv[, 1] * bv[, 2] + K$yy * bv[, 2]^2) /
      pmax(bn^2, .Machine$double.xmin)
    ifelse(bn > 0, sk, pmin(K$xx, K$yy))
  } else {
    eval_spatial(diffusion, cx, cy, name = "diffusion")
  }
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
#' @examples
#' prob <- poisson_problem(pde_mesh_rectangle(nx = 4), source = 1, boundary = dirichlet_bc(0))
#' disc <- discretize(prob)
#' disc
#' dim(disc$A)
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
    D <- supg_diffusion(D, tq, bq)
    if (!is.null(cq)) advt <- list(x = tq * cq * bq$x, y = tq * cq * bq$y)
    gload <- list(x = tq * fq * bq$x, y = tq * fq * bq$y)
  }
  A <- global_matrix(V, local_bilinear(geo, mass = cq, adv = bq, advt = advt, diff = D))
  b <- global_vector(V, local_linear(geo, f = fq, g = gload))
  nat <- natural_bc_terms(V, problem$bcs, t, problem$quad_degree)
  if (!is.null(nat$A)) A <- A + nat$A
  b <- b + nat$b
  dd <- dirichlet_data(V, problem$bcs, t)
  pc <- periodic_constraints(V, problem$bcs)
  structure(list(space = V, A = A, b = b, dirichlet = dd, periodic = pc$P, periodic_col = pc$col,
                 free = setdiff(seq_len(V$ndof), dd$dofs),
                 # Constants are in the kernel exactly when nothing but diffusion,
                 # advection and zero-flux or Neumann data act on u; decided from
                 # the problem data, not from the size of A %*% 1.
                 pure_neumann = !length(dd$dofs) && is.null(cq) &&
                   (is.null(nat$A) || all(nat$A@x == 0)),
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
#' @examples
#' prob <- poisson_problem(pde_mesh_rectangle(nx = 8), source = 1, boundary = dirichlet_bc(0))
#' sol <- solve_pde(prob)
#' summary(sol)
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
  dd <- disc$dirichlet
  free <- disc$free
  A <- disc$A
  bv <- disc$b
  Pm <- disc$periodic
  if (!is.null(Pm)) {
    # Periodic constraints: u = P u_r with one unknown per class of
    # identified degrees of freedom; solve P' A P u_r = P' b.
    col <- disc$periodic_col
    A <- Matrix::drop0(Matrix::crossprod(Pm, A %*% Pm))
    bv <- as.numeric(Matrix::crossprod(Pm, bv))
    rd <- col[dd$dofs]
    keep <- !duplicated(rd)
    dd <- list(dofs = rd[keep], values = dd$values[keep])
    free <- setdiff(seq_len(ncol(Pm)), dd$dofs)
  }
  nd <- nrow(A)
  u <- numeric(nd)
  u[dd$dofs] <- dd$values
  t0 <- proc.time()[["elapsed"]]
  gauge <- NULL
  if (disc$pure_neumann) {
    # Bordered system [A M; M' 0], where column k of M is the mass of the
    # constant on connected component k (one gauge per component). It is
    # nonsingular, and multiplier lambda_k = w_k'b / w_k'm_k vanishes exactly
    # when b is orthogonal to the left null vector w_k of A on component k
    # (w_k = 1 there when A is symmetric, so the test is int f + int g = 0
    # per component; with advection w_k is not constant). lambda_k * sum(m_k)
    # is therefore the compatibility defect of component k.
    one <- assemble_vector(V, 1)
    cc <- mesh_components(V$mesh)
    dc <- dof_components(V, cc)
    if (!is.null(Pm)) {
      if (cc$n > 1L) stopf("Periodic pure Neumann problems need a connected mesh.")
      one <- as.numeric(Matrix::crossprod(Pm, one)); dc <- rep(1L, nd)
    }
    nc <- cc$n
    Mc <- Matrix::sparseMatrix(i = seq_len(nd), j = dc, x = one, dims = c(nd, nc))
    Z <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(nc, nc))
    K <- rbind(cbind(A, Mc), cbind(Matrix::t(Mc), Z))
    s <- solve_linear_system(K, c(bv, numeric(nc)), method = "lu", ...)
    lam <- s$x[nd + seq_len(nc)]
    mass <- as.numeric(Matrix::colSums(Mc))
    comp <- lam * mass
    scale <- as.numeric(tapply(abs(bv), factor(dc, levels = seq_len(nc)), sum))
    rel <- comp / pmax(scale, 1e-300)
    sym <- Matrix::isSymmetric(A, tol = 1e-12)
    # (h / diameter)^2 per component is dimensionless, so the decision does
    # not change with the length unit.
    hrel <- vapply(seq_len(nc), function(k) {
      bb <- apply(V$mesh$nodes[cc$node == k, , drop = FALSE], 2, range)
      mesh_size(V$mesh) / sqrt(sum((bb[2, ] - bb[1, ])^2))
    }, numeric(1))
    tol_c <- if (sym) rep(compatibility_tol, nc) else pmax(compatibility_tol, hrel^2)
    where <- function(k) if (nc > 1L) sprintf(" on connected component %d", k) else ""
    for (kk in seq_len(nc)) {
      if (abs(rel[kk]) <= tol_c[kk]) next
      # With advection the left null vector is only an approximation of the
      # continuous one, so data that are compatible for the PDE leave a
      # defect of the size of the discretization error; only defects above
      # (h / diameter)^2 are reported, and those above 1% are rejected.
      if (sym || abs(rel[kk]) > 0.01 + 10 * tol_c[kk]) {
        stopf(paste0("Pure Neumann problem with incompatible data%s: ",
                     "compatibility defect %s (relative %s). The problem has no solution."),
              where(kk), fmt_num(comp[kk]), fmt_num(rel[kk]))
      }
      warnf(paste0("Pure Neumann problem: the data violate the discrete compatibility condition%s ",
                   "by %s (relative %s); the solution satisfies the equations up to a constant ",
                   "times the multiplier. Check that the data are compatible."),
            where(kk), fmt_num(comp[kk]), fmt_num(rel[kk]))
    }
    u <- s$x[seq_len(nd)]
    gauge <- list(type = "zero mean", components = nc, multiplier = lam,
                  compatibility = comp, relative_compatibility = rel,
                  modified = abs(rel) > tol_c)
    info <- s$info
    # Report the residual of the original equations A u = b; the bordered
    # system is solved to rounding error even when A u = b has no solution.
    info$augmented_relative_residual <- info$relative_residual
    info$residual <- safe_norm(bv - as.numeric(A %*% u))
    nb <- safe_norm(bv)
    info$relative_residual <- if (nb > 0) info$residual / nb else info$residual
    warn_solver(list(status = info$status, relative_residual = info$augmented_relative_residual))
  } else {
    if (!length(free)) {
      info <- list(method = "none", status = "solved", residual = 0, relative_residual = 0,
                   iterations = 0L, n = 0L, nnz = 0L, symmetric = TRUE, tol = NA)
    } else {
      Aff <- A[free, free, drop = FALSE]
      rhs <- bv[free]
      if (length(dd$dofs)) rhs <- rhs - as.numeric(A[free, dd$dofs, drop = FALSE] %*% dd$values)
      dots <- list(...)
      pcs <- dots$preconditioner
      if (inherits(pcs, "pde_preconditioner")) {
        # Tell a mesh-based preconditioner which degrees of freedom are free.
        dots$preconditioner <- structure(list(setup = function(A) pcs$setup(A, free)),
                                         class = "pde_preconditioner")
      }
      s <- do.call(solve_linear_system, c(list(Aff, rhs, method = method), dots))
      u[free] <- s$x
      info <- s$info
      warn_solver(info)
    }
  }
  solve_time <- proc.time()[["elapsed"]] - t0
  diag <- structure(list(
    ndof = V$ndof, free = length(free), nnz = length(A@x),
    assembly_time = disc$assembly_time, solve_time = solve_time,
    solver = info, gauge = gauge,
    full_residual = if (length(free)) sqrt(sum((bv[free] - as.numeric(A[free, , drop = FALSE] %*% u))^2)) else 0
  ), class = "pde_diagnostics")
  if (!is.null(Pm)) u <- as.numeric(Pm %*% u)
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
