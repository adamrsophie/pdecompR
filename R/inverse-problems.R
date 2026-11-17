#' Parameter estimation for PDE models
#'
#' Estimates parameters `theta` of a linear elliptic model from point
#' observations by minimizing
#'
#' `1/2 sum_i ((u_h(x_i; theta) - y_i) / sd_i)^2 + lambda/2 ||theta - theta_ref||^2`
#'
#' subject to the discretized PDE. The gradient is computed with the discrete
#' adjoint method ([pde_sensitivity()]); the optimizer is `stats::optim()`
#' with method `"L-BFGS-B"`, so bounds are supported. Forward solves are
#' cached so the objective and gradient at the same `theta` share one solve.
#'
#' At the estimate, the Gauss-Newton approximation `H = S^T W S + lambda I`
#' (with `S` the forward sensitivities of the predicted observations) gives
#' approximate standard errors. Without regularization they are
#' `sqrt(diag(H^-1))`; with `lambda > 0` they are the propagated measurement
#' noise `sqrt(diag(H^-1 S^T W S H^-1))` (the regularized estimator is biased
#' and this ignores the bias). They reflect measurement noise only, under a
#' local linearization; they do not include discretization error. To expose the latter, the model is re-solved at the
#' estimate on a uniformly refined mesh (when `problem_fun` accepts a mesh)
#' and the change in the predicted observations is reported.
#'
#' A local optimizer is used; for nonconvex problems the result is a local
#' minimum, and no global optimality is claimed.
#'
#' @inheritParams pde_sensitivity
#' @param theta0 Initial parameter values.
#' @param observations A [pde_observations()] object.
#' @param lower,upper Bounds on the parameters.
#' @param regularization Tikhonov weight `lambda`.
#' @param theta_ref Reference parameters of the regularization.
#' @param control Passed to `stats::optim()`.
#' @param derivative Optional exact operator derivatives, see
#'   [pde_sensitivity()] and [affine_derivative()]; indices refer to the full
#'   parameter vector.
#' @return A list of class `pde_estimate` with the estimate, standard errors,
#'   the `fixed` parameters, and `identifiability`: one singular value per
#'   free parameter of the weighted observation sensitivity matrix (zero for
#'   directions the observations cannot see, as when there are fewer
#'   observations than free parameters), its numerical `rank`, the condition
#'   number (infinite when `rank` is below the number of free parameters),
#'   the right singular vectors `directions`, the `null_space` of
#'   unidentified parameter directions, and the correlation matrix (`NULL`
#'   when rank deficient). These describe the information in the data only;
#'   `regularization` can make the optimization well posed but adds no
#'   information.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' pf <- function(theta, mesh) poisson_problem(mesh, diffusion = theta[1], source = 1,
#'                                             boundary = dirichlet_bc(0))
#' truth <- solve_pde(pf(0.5, mesh))
#' xo <- c(0.3, 0.5, 0.7); yo <- c(0.5, 0.5, 0.4)
#' obs <- pde_observations(xo, yo, evaluate(truth$u, xo, yo), sd = 0.01)
#' est <- estimate_parameters(pf, theta0 = 1, observations = obs, mesh = mesh,
#'                            lower = 0.01)
#' est$estimate
#' @export
estimate_parameters <- function(problem_fun, theta0, observations, mesh = NULL,
                                lower = -Inf, upper = Inf, regularization = 0,
                                theta_ref = theta0, control = list(), derivative = NULL) {
  stopifnot(inherits(observations, "pde_observations"))
  theta0 <- as.numeric(theta0)
  np <- length(theta0)
  lower <- rep_len(as.numeric(lower), np); upper <- rep_len(as.numeric(upper), np)
  theta_ref <- rep_len(as.numeric(theta_ref), np)
  if (any(lower > upper)) stopf("Each lower bound must not exceed its upper bound.")
  # Parameters with equal bounds are fixed: optimize and differentiate in the
  # subspace of the free parameters and expand for every model evaluation.
  fixed <- lower == upper
  theta0[fixed] <- lower[fixed]
  ifree <- which(!fixed)
  if (!length(ifree)) stopf("All parameters are fixed by their bounds; there is nothing to estimate.")
  full_theta <- function(th) { out <- theta0_all; out[ifree] <- th; out }
  theta0_all <- theta0
  pf_full <- problem_fun
  problem_fun <- if (length(formals(pf_full)) >= 2L) function(th, m) pf_full(full_theta(th), m)
                 else function(th) pf_full(full_theta(th))
  dfree <- if (is.null(derivative)) NULL else function(th, j) derivative(full_theta(th), ifree[j])
  theta0 <- theta0[ifree]; theta_ref <- theta_ref[ifree]
  lower <- lower[ifree]; upper <- upper[ifree]
  call_pf <- function(th, m = mesh) if (is.null(m)) problem_fun(th) else problem_fun(th, m)
  V <- fe_space(call_pf(theta0)$mesh, call_pf(theta0)$degree)
  obj <- misfit_objective(V, observations)
  cache <- new.env()
  eval_at <- function(theta) {
    if (!is.null(cache$theta) && identical(cache$theta, theta)) return(cache$res)
    s <- pde_sensitivity(problem_fun, theta, mesh = mesh, objective = obj, method = "adjoint",
                         lower = lower, upper = upper, derivative = dfree)
    reg <- regularization * (theta - theta_ref)
    res <- list(value = s$value + 0.5 * sum(reg * (theta - theta_ref)),
                gradient = s$gradient + reg, solution = s$solution)
    cache$theta <- theta; cache$res <- res
    cache$n_solves <- (if (is.null(cache$n_solves)) 0 else cache$n_solves) + 1
    res
  }
  t0 <- proc.time()[["elapsed"]]
  opt <- stats::optim(theta0, function(th) eval_at(th)$value, function(th) eval_at(th)$gradient,
                      method = "L-BFGS-B", lower = lower, upper = upper, control = control)
  est <- opt$par
  fin <- eval_at(est)
  fs <- pde_sensitivity(problem_fun, est, mesh = mesh, method = "forward", lower = lower, upper = upper,
                        derivative = dfree)
  P <- evaluation_matrix(V, observations$x, observations$y)
  So <- as.matrix(P %*% fs$sensitivities)
  W <- 1 / observations$sd^2
  Sw <- So * sqrt(W)
  G <- crossprod(Sw)
  # Identifiability of the free parameters from the weighted observation
  # sensitivities: small singular values (relative to the largest) mark
  # directions the data hardly constrain; correlations near +-1 mark
  # parameters that trade off against each other.
  # The full SVD in parameter space: with fewer observations than free
  # parameters, or dependent sensitivities, the missing singular values are
  # zero and their directions span the null space the data cannot see.
  # Regularization is not information from the data and is not included.
  nfree <- ncol(Sw)
  sv <- svd(Sw, nu = 0, nv = nfree)
  d <- c(sv$d, numeric(nfree - length(sv$d)))
  rank <- sum(d > max(d, 0) * nfree * .Machine$double.eps * 1e3)
  cor_m <- if (rank == nfree) {
    tryCatch({ C <- solve(G); C / sqrt(outer(diag(C), diag(C))) }, error = function(e) NULL)
  }
  ident <- list(singular_values = d, rank = rank, free_parameters = nfree,
                condition = if (rank < nfree) Inf else max(d) / min(d),
                directions = sv$v, null_space = sv$v[, seq_len(nfree) > rank, drop = FALSE],
                correlation = cor_m)
  H <- G + diag(regularization, length(est))
  se <- tryCatch({
    Hi <- solve(H)
    sqrt(diag(if (regularization > 0) Hi %*% G %*% Hi else Hi))
  }, error = function(e) rep(NA_real_, length(est)))
  pred <- as.numeric(P %*% fin$solution$u$values)
  disc <- NULL
  if (length(formals(pf_full)) >= 2L && !is.null(mesh)) {
    fine <- refine_uniform(mesh)
    sf <- solve_pde(call_pf(est, fine))
    pf <- evaluate(sf$u, observations$x, observations$y)
    disc <- list(max_change = max(abs(pf - pred)),
                 relative_to_sd = max(abs(pf - pred) / observations$sd))
  }
  est_full <- full_theta(est); se_full <- rep(NA_real_, np); se_full[ifree] <- se
  grad_full <- rep(NA_real_, np); grad_full[ifree] <- fin$gradient
  structure(list(estimate = est_full, std_error = se_full, fixed = fixed, identifiability = ident,
                 objective = opt$value, gradient = grad_full,
                 convergence = opt$convergence, message = opt$message, counts = opt$counts,
                 forward_solves = cache$n_solves, predicted = pred,
                 residuals = pred - observations$value, solution = fin$solution,
                 discretization_check = disc, time = proc.time()[["elapsed"]] - t0),
            class = "pde_estimate")
}

#' @export
print.pde_estimate <- function(x, ...) {
  cat("PDE parameter estimate (L-BFGS-B, convergence code ", x$convergence, ")\n", sep = "")
  print(data.frame(estimate = x$estimate, std_error = x$std_error), digits = 5)
  cat("  objective ", fmt_num(x$objective), ", gradient norm ", fmt_num(sqrt(sum(x$gradient^2))),
      ", forward solves ", x$forward_solves, "\n", sep = "")
  cat("  standard errors: measurement noise only (Gauss-Newton linearization)\n")
  if (!is.null(x$identifiability)) {
    id <- x$identifiability
    cat("  identifiability: rank ", id$rank, " of ", id$free_parameters, " free parameters, condition number ",
        fmt_num(id$condition, 3), "\n", sep = "")
  }
  if (!is.null(x$discretization_check)) {
    cat("  discretization check: predictions change by up to ", fmt_num(x$discretization_check$max_change),
        " (", fmt_num(x$discretization_check$relative_to_sd), " sd) on a refined mesh\n", sep = "")
  }
  invisible(x)
}

#' Distributed optimal control
#'
#' Solves the linear-quadratic optimal control problem
#'
#' `min_q 1/2 ||u - u_target||^2 + alpha/2 ||q||^2`
#'
#' subject to `-div(k grad u) = f + q` with the given boundary conditions,
#' where the control `q` is a P1 function and norms are `L2` norms (computed
#' with the mass matrix). The target enters through its nodal interpolant, so
#' the reported tracking term is `||u - I_h u_target||`; features of the
#' target smaller than the mesh size are not seen. The problem is strictly convex for `alpha > 0`.
#'
#' Without bounds the first-order optimality (KKT) system is solved directly
#' with sparse LU. With bounds `lower <= q <= upper`, the reduced problem is
#' solved with `stats::optim()` (`"L-BFGS-B"`) using adjoint gradients and a
#' single Cholesky factorization of the state operator.
#'
#' @param mesh A [pde_mesh()].
#' @param target Desired state `u_target(x, y)`.
#' @param alpha Control cost (positive).
#' @param diffusion Diffusion coefficient.
#' @param source Fixed source `f`.
#' @param bcs Boundary conditions of the state equation.
#' @param lower,upper Optional control bounds (numbers).
#' @param control Passed to `stats::optim()`.
#' @return A list of class `pde_control` with the `control` and `state`
#'   functions, the `objective`, the gradient norm and solver details.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' oc <- optimal_control(mesh, target = function(x, y) sin(pi * x) * sin(pi * y),
#'                       alpha = 1e-4)
#' oc
#' @export
optimal_control <- function(mesh, target, alpha, diffusion = 1, source = 0,
                            bcs = dirichlet_bc(0), lower = NULL, upper = NULL,
                            control = list()) {
  check_number(alpha, "alpha", positive = TRUE)
  disc <- discretize(elliptic_problem(mesh, diffusion = diffusion, source = source, bcs = bcs))
  V <- disc$space
  if (disc$pure_neumann) stopf("The state equation must have a unique solution.")
  M <- assemble_matrix(V, "mass")
  ud <- interpolate(target, V)$values
  # Constrained (periodic) state unknowns u = P u_r; the control stays a
  # full finite-element function.
  cs <- constrained_system(disc); Pm <- cs$P
  f <- cs$free; d <- cs$dirichlet$dofs; uD <- cs$dirichlet$values
  A <- cs$A
  Aff <- A[f, f, drop = FALSE]
  r0 <- cs$b[f]
  if (length(d)) r0 <- r0 - as.numeric(A[f, d, drop = FALSE] %*% uD)
  PM <- Matrix::crossprod(Pm, M)          # P' M
  Mr <- PM %*% Pm                         # P' M P
  Mf <- PM[f, , drop = FALSE]
  bad <- new.env(); bad$n <- 0L
  note <- function(info) if (!solve_ok(info)) { bad$n <- bad$n + 1L; bad$info <- info }
  # Paired state cache: the objective and gradient at the same control share
  # one state solve.
  last <- new.env()
  state <- function(q, slv) {
    if (!is.null(last$q) && identical(last$q, q)) return(last$u)
    ur <- numeric(ncol(Pm)); ur[d] <- uD
    sv <- slv$solve(r0 + as.numeric(Mf %*% q))
    note(sv$info)
    ur[f] <- sv$x
    u <- as.numeric(Pm %*% ur)
    last$q <- q; last$u <- u
    u
  }
  J <- function(u, q) {
    e <- u - ud
    0.5 * sum(e * as.numeric(M %*% e)) + 0.5 * alpha * sum(q * as.numeric(M %*% q))
  }
  t0 <- proc.time()[["elapsed"]]
  slv <- linear_solver(Aff)
  # The adjoint solves with A^T (A is symmetric for the supported scalar
  # symmetric diffusion; the transpose keeps the gradient exact otherwise).
  slv_t <- if (isTRUE(slv$info$symmetric)) slv else linear_solver(Matrix::t(Aff))
  grad_fun <- function(q) {
    u <- state(q, slv)
    av <- slv_t$solve(as.numeric(PM %*% (u - ud))[f]); note(av$info)
    lam <- av$x
    as.numeric(Matrix::crossprod(Mf, lam)) + alpha * as.numeric(M %*% q)
  }
  bounded <- !is.null(lower) || !is.null(upper)
  if (!bounded) {
    nf <- length(f); n <- V$ndof
    Z <- function(r, c) Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(r, c))
    Mff <- Mr[f, f, drop = FALSE]
    rhs1 <- as.numeric(PM %*% ud)[f]
    if (length(d)) rhs1 <- rhs1 - as.numeric(Mr[f, d, drop = FALSE] %*% uD)
    KKT <- rbind(cbind(Mff, Z(nf, n), -Matrix::t(Aff)),
                 cbind(Z(n, nf), alpha * M, Matrix::t(Mf)),
                 cbind(Aff, -Mf, Z(nf, nf)))
    s <- solve_linear_system(KKT, c(rhs1, numeric(n), r0), method = "lu")
    note(s$info)
    q <- s$x[nf + seq_len(n)]
    method <- "KKT (sparse LU)"; iters <- NA; conv <- 0L
  } else {
    lo <- if (is.null(lower)) -Inf else lower; hi <- if (is.null(upper)) Inf else upper
    q0 <- pmin(pmax(numeric(V$ndof), lo), hi)
    opt <- stats::optim(q0, function(q) J(state(q, slv), q), grad_fun, method = "L-BFGS-B",
                        lower = lo, upper = hi, control = utils::modifyList(list(maxit = 500), control))
    q <- opt$par
    method <- "reduced L-BFGS-B (adjoint gradient)"; iters <- opt$counts; conv <- opt$convergence
  }
  u <- state(q, slv)
  g <- grad_fun(q)
  if (bounded) {
    lo <- if (is.null(lower)) -Inf else lower; hi <- if (is.null(upper)) Inf else upper
    act <- (q <= lo + 1e-10 & g > 0) | (q >= hi - 1e-10 & g < 0)
    g[act] <- 0
  }
  if (bad$n) warn_solver(bad$info, sprintf("state, adjoint or KKT solve (%d time(s))", bad$n))
  structure(list(control = fe_function(V, q, "q"), state = fe_function(V, u, "u"),
                 target = fe_function(V, ud, "u_target"), objective = J(u, q),
                 gradient_norm = sqrt(sum(g^2)), method = method, iterations = iters,
                 convergence = conv, solver_failures = bad$n, time = proc.time()[["elapsed"]] - t0),
            class = "pde_control")
}

#' @export
print.pde_control <- function(x, ...) {
  cat("Distributed optimal control (", x$method, ")\n", sep = "")
  cat("  objective ", fmt_num(x$objective, 6), ", projected gradient norm ",
      fmt_num(x$gradient_norm, 3), "\n", sep = "")
  cat("  control range [", fmt_num(min(x$control$values)), ", ", fmt_num(max(x$control$values)),
      "]\n", sep = "")
  invisible(x)
}
