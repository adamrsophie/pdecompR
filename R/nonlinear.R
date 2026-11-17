#' Nonlinear stationary problems
#'
#' `nonlinear_problem()` defines the quasilinear problem
#'
#' `-div(k(u) grad u) + g(u) = f`
#'
#' and `p_laplace_problem()` the regularized p-Laplacian
#'
#' `-div((|grad u|^2 + eps^2)^((p - 2) / 2) grad u) = f`.
#'
#' Both are solved by Newton's method on the residual `F(u) = 0` of the weak
#' form, with the exact (consistent) Jacobian assembled from the derivatives
#' `k'(u)` and `g'(u)` (finite differences if not supplied), and a
#' backtracking line search on the residual norm. The convergence history is
#' returned; Newton is not assumed to converge, and failure is reported.
#'
#' The initial guess is the solution of the linear problem with `k(0)` and
#' `g'(0)` (for the p-Laplacian, the Laplace problem), unless `initial` is
#' given.
#'
#' @inheritParams elliptic_problem
#' @param diffusion `k(u)`: a number or `function(u)`.
#' @param diffusion_derivative `k'(u)` as `function(u)`, or `NULL`.
#' @param reaction `g(u)`: `NULL` or `function(u)`.
#' @param reaction_derivative `g'(u)` as `function(u)`, or `NULL`.
#' @param p Exponent of the p-Laplacian (`p > 1`).
#' @param epsilon Regularization parameter of the p-Laplacian.
#' @param initial Optional initial guess (number, function or [fe_function()]).
#' @return An object of class `nonlinear_problem`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' prob <- nonlinear_problem(mesh, diffusion = function(u) 1 + u^2,
#'                           diffusion_derivative = function(u) 2 * u,
#'                           source = 10, bcs = dirichlet_bc(0))
#' sol <- solve_pde(prob)
#' sol$newton
#' @references Deuflhard, P. (2004). *Newton Methods for Nonlinear Problems*.
#'   Springer.
#' @export
nonlinear_problem <- function(mesh, diffusion = 1, diffusion_derivative = NULL,
                              reaction = NULL, reaction_derivative = NULL, source = 0,
                              bcs = list(), degree = 1, initial = NULL, quad_degree = NULL) {
  stopifnot(inherits(mesh, "pde_mesh"))
  bcs <- normalize_bcs(bcs)
  no_periodic(bcs, "nonlinear problems")
  check_bc_markers(bcs, mesh)
  if (!(is.function(diffusion) || (is.numeric(diffusion) && length(diffusion) == 1L))) {
    stopf("`diffusion` must be a number or a function of u.")
  }
  if (!is.null(reaction) && !is.function(reaction)) stopf("`reaction` must be a function of u.")
  check_coefficient(source, "source")
  structure(list(mesh = mesh, type = "quasilinear", diffusion = diffusion,
                 diffusion_derivative = diffusion_derivative, reaction = reaction,
                 reaction_derivative = reaction_derivative, source = source, bcs = bcs,
                 degree = check_count(degree, "degree"), initial = initial,
                 quad_degree = quad_degree),
            class = c("nonlinear_problem", "pde_problem"))
}

#' @rdname nonlinear_problem
#' @export
p_laplace_problem <- function(mesh, p, source = 0, bcs = list(), epsilon = 1e-6, degree = 1,
                              initial = NULL, quad_degree = NULL) {
  check_number(p, "p")
  if (p <= 1) stopf("`p` must be larger than 1.")
  check_number(epsilon, "epsilon", positive = TRUE)
  prob <- nonlinear_problem(mesh, source = source, bcs = bcs, degree = degree,
                            initial = initial, quad_degree = quad_degree)
  prob$type <- "p-laplace"; prob$p <- p; prob$epsilon <- epsilon
  prob
}

fd_deriv <- function(f) function(u) {
  h <- 1e-6 * pmax(1, abs(u))
  (f(u + h) - f(u - h)) / (2 * h)
}

as_fun_u <- function(k) if (is.function(k)) k else function(u) rep(k, length(u))

# Evaluate a function of u at quadrature values, recycling constant results.
eval_u <- function(f, U, name) {
  v <- f(as.vector(U))
  if (length(v) == 1L) v <- rep(v, length(U))
  if (length(v) != length(U) || !is.numeric(v)) stopf("`%s` must return one value per input value.", name)
  if (!all(is.finite(v))) stop_nonfinite("`%s` returned non-finite values.", name)
  matrix(v, nrow(U))
}

#' @rdname solve_pde
#' @param tol,maxit Newton tolerance on the free residual norm (relative to the
#'   initial residual, with an absolute floor of `1e-12`) and maximum number
#'   of Newton iterations.
#' @export
solve_pde.nonlinear_problem <- function(problem, tol = 1e-10, maxit = 30, method = "auto", ...) {
  t0 <- proc.time()[["elapsed"]]
  check_tol(tol, "tol"); maxit <- check_count(maxit, "maxit")
  V <- fe_space(problem$mesh, problem$degree)
  geo <- element_geometry(V, problem$quad_degree)
  fq <- coef_qp(problem$source, geo, name = "source")
  nat <- natural_bc_terms(V, problem$bcs)
  b_nat <- global_vector(V, local_linear(geo, f = fq)) + nat$b
  dd <- dirichlet_data(V, problem$bcs)
  free <- setdiff(seq_len(V$ndof), dd$dofs)
  plap <- identical(problem$type, "p-laplace")
  k <- as_fun_u(problem$diffusion)
  dk <- if (!is.null(problem$diffusion_derivative)) problem$diffusion_derivative
        else if (is.function(problem$diffusion)) fd_deriv(problem$diffusion) else function(u) 0 * u
  g <- problem$reaction
  dg <- if (is.null(g)) NULL else if (!is.null(problem$reaction_derivative)) problem$reaction_derivative else fd_deriv(g)
  # Initial guess
  if (!is.null(problem$initial)) {
    u <- initial_values(problem$initial, V)
  } else {
    # Linearize at the mean Dirichlet value; fall back to unit diffusion when
    # k is not positive there (for example k(u) = u^2 at u = 0).
    u0 <- if (length(dd$values)) mean(dd$values) else 0
    k0 <- if (plap) 1 else k(u0)
    if (!is.finite(k0) || k0 <= 0) k0 <- 1
    lin <- elliptic_problem(problem$mesh, diffusion = k0,
                            reaction = if (is.null(dg)) 0 else max(0, dg(0)),
                            source = if (is.null(g)) problem$source else {
                              s0 <- problem$source; g0 <- g(0)
                              function(x, y) eval_spatial(s0, x, y) - g0
                            },
                            bcs = problem$bcs, degree = problem$degree,
                            quad_degree = problem$quad_degree)
    u <- solve_pde(lin)$u$values
  }
  u[dd$dofs] <- dd$values
  if (!length(dd$dofs)) {
    # Without Dirichlet data the solution is unique only through a reaction or
    # Robin term. A reaction that vanishes with its derivative at the initial
    # state, or an all-zero Robin matrix, counts as absent (a narrow check, not
    # a general uniqueness proof).
    robin <- !is.null(nat$A) && any(nat$A@x != 0)
    react <- FALSE
    if (!is.null(g)) {
      uq <- fe_qp(u, geo)$value
      react <- any(eval_u(g, uq, "reaction") != 0) || any(eval_u(dg, uq, "reaction_derivative") != 0)
    }
    if (!robin && !react) {
      stopf(paste0("The nonlinear problem has no Dirichlet, Robin or reaction term, so its solution ",
                   "is not unique (pure Neumann); add a Dirichlet condition."))
    }
  }
  residual_jac <- function(u, jacobian = TRUE) {
    q <- fe_qp(u, geo)
    if (plap) {
      s <- q$dx^2 + q$dy^2 + problem$epsilon^2
      kq <- s^((problem$p - 2) / 2)
    } else {
      kq <- eval_u(k, q$value, "diffusion")
    }
    flux <- list(x = kq * q$dx, y = kq * q$dy)
    gq <- if (is.null(g)) NULL else eval_u(g, q$value, "reaction")
    r <- global_vector(V, local_linear(geo, f = gq, g = flux)) - b_nat
    if (!is.null(nat$A)) r <- r + as.numeric(nat$A %*% u)
    if (!jacobian) return(list(r = r))
    if (plap) {
      c2 <- (problem$p - 2) * s^((problem$p - 4) / 2)
      D <- list(xx = kq + c2 * q$dx * q$dx, xy = c2 * q$dx * q$dy,
                yx = c2 * q$dy * q$dx, yy = kq + c2 * q$dy * q$dy)
      A <- local_bilinear(geo, diff = D)
    } else {
      dkq <- eval_u(dk, q$value, "diffusion_derivative")
      mass <- if (is.null(dg)) NULL else eval_u(dg, q$value, "reaction_derivative")
      A <- local_bilinear(geo, mass = mass, advt = list(x = dkq * q$dx, y = dkq * q$dy), diff = kq)
    }
    J <- global_matrix(V, A)
    if (!is.null(nat$A)) J <- J + nat$A
    list(r = r, J = J)
  }
  hist <- data.frame(iteration = integer(), residual = numeric(), step = numeric())
  rj <- residual_jac(u)
  r0 <- sqrt(sum(rj$r[free]^2)); rn <- r0
  hist[1, ] <- list(0L, r0, NA)
  it <- 0L
  converged <- rn <= 1e-12
  status <- "not converged"
  while (!converged && it < maxit) {
    it <- it + 1L
    ls <- solve_linear_system(rj$J[free, free, drop = FALSE], -rj$r[free], method = method)
    du <- ls$x
    if (!ls$info$status %in% c("solved", "converged") || !all(is.finite(du))) {
      status <- "linear solve failed"; break
    }
    lam <- 1; accepted <- FALSE
    while (lam >= 1e-6) {
      ut <- u; ut[free] <- u[free] + lam * du
      rt <- trial_eval(function() residual_jac(ut, jacobian = FALSE)$r)
      rtn <- if (is.null(rt)) Inf else sqrt(sum(rt[free]^2))
      if (is.finite(rtn) && rtn < (1 - 1e-4 * lam) * rn) { accepted <- TRUE; break }
      lam <- lam / 2
    }
    if (!accepted) { status <- "line search failed"; break }
    u <- ut
    rj <- residual_jac(u)
    rn <- sqrt(sum(rj$r[free]^2))
    hist[it + 1L, ] <- list(it, rn, lam)
    converged <- rn <= max(tol * r0, 1e-12)
  }
  if (converged) status <- "converged"
  if (!converged) {
    warnf("Newton stopped (%s) after %d iterations (residual %s, initial %s).",
          status, it, fmt_num(rn), fmt_num(r0))
  }
  diag <- structure(list(ndof = V$ndof, free = length(free), nnz = length(rj$J@x),
                         assembly_time = NA_real_, solve_time = proc.time()[["elapsed"]] - t0,
                         solver = list(method = "newton", status = status,
                                       residual = rn, relative_residual = rn / max(r0, 1e-300),
                                       iterations = it),
                         gauge = NULL),
                    class = "pde_diagnostics")
  structure(list(u = fe_function(V, u), space = V, diagnostics = diag, newton = hist,
                 converged = converged, problem = problem), class = "pde_solution")
}
