#' Parameter sensitivities of linear elliptic problems
#'
#' For a parameterized linear problem, the reduced discrete system after
#' Dirichlet elimination is `K(theta) u_f = r(theta)` with Dirichlet values
#' `u_d(theta)`. Differentiating gives the forward sensitivity equation
#'
#' `K du_f/dtheta_j = dr/dtheta_j - dK/dtheta_j u_f`,
#'
#' valid when `K` is invertible and the data are differentiable. For a scalar
#' objective `J(u)` the discrete adjoint method solves `K^T lambda = dJ/du_f`
#' once and returns `dJ/dtheta_j = dJ/du_d . du_d/dtheta_j +
#' lambda . (dr/dtheta_j - dK/dtheta_j u_f)`, at a cost independent of the
#' number of parameters (apart from the operator derivatives).
#'
#' These are discrete sensitivities (differentiate, then the discrete system
#' is exact); they coincide with continuous adjoints only up to
#' discretization error.
#'
#' The operator derivatives `dK/dtheta_j`, `dr/dtheta_j` and `du_d/dtheta_j`
#' are obtained by central differences of the assembled operators with step
#' `h_j = h max(1, |theta_j|)`. This is exact up to rounding when the
#' coefficients are affine in `theta`, and second-order accurate otherwise.
#' The solves themselves are not differenced. Supply `derivative` for exact
#' operator derivatives (for example of an affine family
#' `A(theta) = A0 + sum(theta_j A_j)`, where `derivative(theta, j)` returns
#' `A_j`); check either route with [taylor_test()]. Each parameter's
#' derivative is formed, used and discarded in turn, so memory does not grow
#' with the number of parameters.
#'
#' @param problem_fun `function(theta)` or `function(theta, mesh)` returning
#'   an [elliptic_problem()].
#' @param theta Parameter vector.
#' @param mesh Optional mesh passed as second argument of `problem_fun`.
#' @param objective Optional `function(u)` of the full degree-of-freedom
#'   vector returning `list(value = , gradient = )`, for example from
#'   [misfit_objective()].
#' @param method `"forward"` (all solution sensitivities) or `"adjoint"`
#'   (objective gradient only).
#' @param h Relative difference step for the operator derivatives.
#' @param lower,upper Parameter bounds. Difference steps are shortened or
#'   made one-sided so that `problem_fun` is never called outside them.
#' @param derivative Optional exact derivative `function(theta, j)` returning
#'   `list(A = , b = , dirichlet = )`: the derivatives of the full system
#'   matrix and load vector (before Dirichlet elimination, as in
#'   [discretize()]) and, optionally, of the Dirichlet values with respect to
#'   `theta[j]`.
#' @return A list of class `pde_sensitivity` with the `solution`, the
#'   `sensitivities` matrix (forward), and the objective `value` and
#'   `gradient` when `objective` is given.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' pf <- function(theta) poisson_problem(mesh, diffusion = theta[1], source = 1,
#'                                       boundary = dirichlet_bc(0))
#' s <- pde_sensitivity(pf, theta = 2)
#' # exact: u = u1 / theta, so du/dtheta = -u / theta
#' max(abs(s$sensitivities[, 1] + s$solution$u$values / 2))
#' @references Hinze, M., Pinnau, R., Ulbrich, M. and Ulbrich, S. (2009).
#'   *Optimization with PDE Constraints*. Springer.
#' @export
pde_sensitivity <- function(problem_fun, theta, mesh = NULL, objective = NULL,
                            method = c("forward", "adjoint"), h = 1e-6, lower = -Inf, upper = Inf,
                            derivative = NULL) {
  method <- match.arg(method)
  if (method == "adjoint" && is.null(objective)) stopf("The adjoint method needs an `objective`.")
  if (!is.null(derivative) && !is.function(derivative)) stopf("`derivative` must be a function(theta, j).")
  check_tol(h, "h")
  theta <- as.numeric(theta)
  if (!length(theta) || any(!is.finite(theta))) stopf("`theta` must be a nonempty finite vector.")
  call_pf <- function(th) if (is.null(mesh)) problem_fun(th) else problem_fun(th, mesh)
  base <- reduced_system(call_pf(theta))
  u <- base$u
  p <- length(theta)
  # Work on the constrained (periodic) unknowns; P lifts them to all dofs.
  Pm <- base$P
  red <- function(g) as.numeric(Matrix::crossprod(Pm, g))
  uf <- base$ur[base$free]
  ud_full <- base$ud
  # Derivatives of the reduced operators for parameter j, computed when
  # needed and discarded after use, so at most one derivative matrix is held
  # at a time.
  deriv_j <- function(j) {
    if (!is.null(derivative)) {
      dj <- derivative(theta, j)
      dA <- dj$A; db <- dj$b
      n <- length(u)
      if (!inherits(dA, c("Matrix", "matrix")) || !all(dim(dA) == n) || length(db) != n) {
        stopf("`derivative(theta, j)` must return list(A = n by n matrix, b = length-n vector).")
      }
      dA <- Matrix::crossprod(Pm, dA %*% Pm); db <- red(db)
      dud <- if (is.null(dj$dirichlet)) numeric(length(base$dofs)) else as.numeric(dj$dirichlet)
      # A derivative of the full Dirichlet data (as affine_derivative()
      # returns) is reduced by the same periodic classes as the data; values
      # within one class must agree.
      if (!is.null(base$col) && length(dud) == length(base$full_dofs) && length(dud) != length(base$dofs)) {
        dud <- reduce_prescribed(list(dofs = base$full_dofs, values = dud), base$col)$values
      }
      if (length(dud) != length(base$dofs)) stopf("`derivative(theta, j)$dirichlet` must have one value per Dirichlet dof.")
      f <- base$free; d <- base$dofs
      dK <- dA[f, f, drop = FALSE]
      dr <- as.numeric(db[f])
      if (length(d)) {
        dr <- dr - as.numeric(dA[f, d, drop = FALSE] %*% ud_full) -
          as.numeric(base$A_fd %*% dud)
      }
      return(list(dK = dK, dr = dr, dud = dud))
    }
    hj <- h * max(1, abs(theta[j]))
    lo <- rep_len(lower, p)[j]; hi <- rep_len(upper, p)[j]
    # Shrink the step or switch to a one-sided difference to stay in bounds.
    up <- min(hj, hi - theta[j]); dn <- min(hj, theta[j] - lo)
    if (up <= 0 && dn <= 0) stopf("Parameter %d has no room between its bounds for a difference step.", j)
    if (up > 0 && dn > 0) { up <- dn <- min(up, dn) } else if (up <= 0) up <- 0 else dn <- 0
    tp <- theta; tp[j] <- tp[j] + up
    tm <- theta; tm[j] <- tm[j] - dn
    hj <- (up + dn) / 2
    sp <- reduced_system(call_pf(tp), solve = FALSE)
    sm <- reduced_system(call_pf(tm), solve = FALSE)
    if (!identical(sp$dofs, base$dofs) || !identical(sm$dofs, base$dofs)) {
      stopf("The Dirichlet degrees of freedom must not depend on theta.")
    }
    list(dK = (sp$K - sm$K) / (2 * hj), dr = (sp$r - sm$r) / (2 * hj), dud = (sp$ud - sm$ud) / (2 * hj))
  }
  out <- list(solution = base$solution, theta = theta, method = method,
              derivatives = if (is.null(derivative)) "finite difference" else "exact")
  if (method == "forward") {
    S <- matrix(0, ncol(Pm), p)
    for (j in seq_len(p)) {
      dj <- deriv_j(j)
      sj <- base$solver$solve(dj$dr - as.numeric(dj$dK %*% uf))
      warn_solver(sj$info, "sensitivity solve")
      S[base$free, j] <- sj$x
      if (length(base$dofs)) S[base$dofs, j] <- dj$dud
    }
    S <- as.matrix(Pm %*% S)
    out$sensitivities <- S
    if (!is.null(objective)) {
      ob <- objective(u)
      out$value <- ob$value
      out$gradient <- as.numeric(crossprod(S, ob$gradient))
    }
  } else {
    ob <- objective(u)
    gr_r <- red(ob$gradient)
    KT <- if (isTRUE(base$solver$info$symmetric)) base$solver else linear_solver(Matrix::t(base$K))
    as <- KT$solve(gr_r[base$free])
    warn_solver(as$info, "adjoint solve")
    lam <- as$x
    gr <- vapply(seq_len(p), function(j) {
      dj <- deriv_j(j)
      sum(lam * (dj$dr - as.numeric(dj$dK %*% uf))) +
        if (length(base$dofs)) sum(gr_r[base$dofs] * dj$dud) else 0
    }, numeric(1))
    out$value <- ob$value
    out$gradient <- gr
    out$adjoint <- lam
  }
  structure(out, class = "pde_sensitivity")
}

#' Exact derivatives of an affine problem family
#'
#' For problems whose assembled system is affine in the parameters,
#' `A(theta) = A0 + sum(theta_j A_j)` and `b(theta) = b0 + sum(theta_j b_j)`
#' (for example a diffusion, reaction or source that is a linear
#' combination of fixed fields), the components are the unit differences
#' `A_j = A(theta + e_j) - A(theta)`. `affine_derivative()` assembles them once
#' and returns a `derivative` function for [pde_sensitivity()] and
#' [taylor_test()]; later calls cost no assembly. The result is exact only
#' for affine families: verify with [taylor_test()].
#'
#' @param problem_fun,mesh As in [pde_sensitivity()].
#' @param theta Parameter vector at which the components are assembled.
#' @return A `function(theta, j)` returning `list(A = , b = , dirichlet = )`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 6)
#' pf <- function(theta) poisson_problem(mesh, diffusion = theta[1], source = theta[2],
#'                                       boundary = dirichlet_bc(0))
#' dv <- affine_derivative(pf, c(1, 1))
#' obs <- pde_observations(0.5, 0.5, 0.1)
#' pde_sensitivity(pf, c(2, 3), objective = misfit_objective(fe_space(mesh), obs),
#'                 method = "adjoint", derivative = dv)$gradient
#' @export
affine_derivative <- function(problem_fun, theta, mesh = NULL) {
  theta <- as.numeric(theta)
  call_pf <- function(th) if (is.null(mesh)) problem_fun(th) else problem_fun(th, mesh)
  d0 <- discretize(call_pf(theta))
  comps <- lapply(seq_along(theta), function(j) {
    tj <- theta; tj[j] <- tj[j] + 1
    dj <- discretize(call_pf(tj))
    if (!identical(dj$dirichlet$dofs, d0$dirichlet$dofs)) {
      stopf("The Dirichlet degrees of freedom must not depend on theta.")
    }
    list(A = Matrix::drop0(dj$A - d0$A), b = dj$b - d0$b,
         dirichlet = dj$dirichlet$values - d0$dirichlet$values)
  })
  function(theta, j) comps[[j]]
}

#' Taylor remainder test of an objective gradient
#'
#' Checks a gradient `g` of `J(theta)` by the remainders
#' `|J(theta + e d) - J(theta) - e g . d|` for decreasing `e`. They must
#' decrease at second order (observed rates near 2) for a correct gradient
#' and only at first order for a wrong one.
#'
#' @param problem_fun,mesh,objective,lower,upper,derivative As in
#'   [pde_sensitivity()].
#' @param theta Parameter vector.
#' @param direction Perturbation direction (default: all ones).
#' @param eps Step sizes, decreasing.
#' @param method Gradient method, `"adjoint"` or `"forward"`.
#' @return A data frame with `eps`, `remainder` and observed `rate`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 6)
#' pf <- function(theta) poisson_problem(mesh, diffusion = theta[1], source = theta[2],
#'                                       boundary = dirichlet_bc(0))
#' V <- fe_space(mesh)
#' obs <- pde_observations(0.5, 0.5, 0.1)
#' taylor_test(pf, c(1, 1), objective = misfit_objective(V, obs))
#' @export
taylor_test <- function(problem_fun, theta, objective, direction = NULL, mesh = NULL,
                        eps = 10^-(1:4), method = c("adjoint", "forward"), derivative = NULL,
                        lower = -Inf, upper = Inf) {
  method <- match.arg(method)
  theta <- as.numeric(theta)
  d <- if (is.null(direction)) rep(1, length(theta)) else as.numeric(direction)
  if (length(d) != length(theta)) stopf("`direction` must have one entry per parameter.")
  s0 <- pde_sensitivity(problem_fun, theta, mesh = mesh, objective = objective, method = method,
                        derivative = derivative, lower = lower, upper = upper)
  call_pf <- function(th) if (is.null(mesh)) problem_fun(th) else problem_fun(th, mesh)
  Jt <- function(th) objective(solve_pde(call_pf(th))$u$values)$value
  rem <- vapply(eps, function(e) abs(Jt(theta + e * d) - s0$value - e * sum(s0$gradient * d)), numeric(1))
  data.frame(eps = eps, remainder = rem,
             rate = c(NA, log(rem[-length(rem)] / rem[-1]) / log(eps[-length(eps)] / eps[-1])))
}

#' @export
print.pde_sensitivity <- function(x, ...) {
  cat("<pde_sensitivity>", x$method, "sensitivities for", length(x$theta), "parameter(s)\n")
  if (!is.null(x$gradient)) cat("  objective", fmt_num(x$value), " gradient", fmt_num(x$gradient), "\n")
  invisible(x)
}

reduced_system <- function(problem, solve = TRUE) {
  if (!inherits(problem, "elliptic_problem")) stopf("`problem_fun` must return an elliptic_problem.")
  if (problem$stabilization != "none") stopf("Sensitivities with stabilization are not supported.")
  disc <- discretize(problem)
  if (disc$pure_neumann) stopf("Sensitivities of pure Neumann problems are not supported.")
  cs <- constrained_system(disc)
  f <- cs$free; d <- cs$dirichlet$dofs; A <- cs$A
  K <- A[f, f, drop = FALSE]
  r <- cs$b[f]
  if (length(d)) r <- r - as.numeric(A[f, d, drop = FALSE] %*% cs$dirichlet$values)
  out <- list(K = K, r = r, ud = cs$dirichlet$values, free = f, dofs = d,
              A_fd = A[f, d, drop = FALSE], P = cs$P, full_dofs = disc$dirichlet$dofs, col = disc$periodic_col)
  if (solve) {
    slv <- linear_solver(K)
    ur <- numeric(ncol(cs$P))
    ur[d] <- cs$dirichlet$values
    ts <- proc.time()[["elapsed"]]
    sv <- slv$solve(r)
    warn_solver(sv$info, "state solve")
    ur[f] <- sv$x
    u <- as.numeric(cs$P %*% ur)
    out$u <- u; out$ur <- ur
    out$solver <- slv
    diag <- structure(list(ndof = disc$space$ndof, free = length(f), nnz = length(disc$A@x),
                           solver = sv$info, assembly_time = disc$assembly_time,
                           solve_time = proc.time()[["elapsed"]] - ts),
                      class = "pde_diagnostics")
    out$solution <- structure(list(u = fe_function(disc$space, u), space = disc$space,
                                   problem = problem, diagnostics = diag), class = "pde_solution")
  }
  out
}

#' Point observations and data misfit
#'
#' `pde_observations()` stores point measurements `y_i` at `(x_i, y_i)` with
#' standard deviations `sd_i`. `misfit_objective()` returns the weighted
#' least-squares misfit `J(u) = 1/2 sum_i ((P u)_i - y_i)^2 / sd_i^2` and its
#' gradient `P^T W (P u - y)`, where `P` is the [evaluation_matrix()].
#'
#' @param x,y Observation coordinates.
#' @param values Observed values.
#' @param sd Measurement standard deviations (recycled).
#' @param space An [fe_space()].
#' @param observations A `pde_observations` object.
#' @return `pde_observations()` returns a data frame of class
#'   `pde_observations`; `misfit_objective()` returns a function of the
#'   degree-of-freedom vector.
#' @examples
#' obs <- pde_observations(c(0.3, 0.7), c(0.5, 0.5), c(0.1, 0.2), sd = 0.01)
#' V <- fe_space(pde_mesh_rectangle(nx = 4))
#' J <- misfit_objective(V, obs)
#' J(interpolate(function(x, y) x, V)$values)$value
#' @export
pde_observations <- function(x, y, values, sd = 1) {
  n <- length(x)
  if (length(y) != n || length(values) != n) stopf("`x`, `y` and `values` must have equal length.")
  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(values) ||
      any(!is.finite(c(x, y, values)))) {
    stopf("`x`, `y` and `values` must be finite numbers.")
  }
  sd <- rep_len(sd, n)
  if (any(!is.finite(sd) | sd <= 0)) stopf("`sd` must be positive.")
  structure(data.frame(x = x, y = y, value = values, sd = sd), class = c("pde_observations", "data.frame"))
}

#' @rdname pde_observations
#' @export
misfit_objective <- function(space, observations) {
  P <- evaluation_matrix(space, observations$x, observations$y)
  w <- 1 / observations$sd^2
  y <- observations$value
  function(u) {
    r <- as.numeric(P %*% u) - y
    list(value = 0.5 * sum(w * r^2), gradient = as.numeric(Matrix::crossprod(P, w * r)))
  }
}
