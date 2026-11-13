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
#' The solves themselves are not differenced.
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
#' @export
pde_sensitivity <- function(problem_fun, theta, mesh = NULL, objective = NULL,
                            method = c("forward", "adjoint"), h = 1e-6) {
  method <- match.arg(method)
  if (method == "adjoint" && is.null(objective)) stopf("The adjoint method needs an `objective`.")
  theta <- as.numeric(theta)
  call_pf <- function(th) if (is.null(mesh)) problem_fun(th) else problem_fun(th, mesh)
  base <- reduced_system(call_pf(theta))
  u <- base$u
  p <- length(theta)
  dK <- vector("list", p); dr <- vector("list", p); dud <- vector("list", p)
  for (j in seq_len(p)) {
    hj <- h * max(1, abs(theta[j]))
    tp <- theta; tp[j] <- tp[j] + hj
    tm <- theta; tm[j] <- tm[j] - hj
    sp <- reduced_system(call_pf(tp), solve = FALSE)
    sm <- reduced_system(call_pf(tm), solve = FALSE)
    if (!identical(sp$dofs, base$dofs) || !identical(sm$dofs, base$dofs)) {
      stopf("The Dirichlet degrees of freedom must not depend on theta.")
    }
    dK[[j]] <- (sp$K - sm$K) / (2 * hj)
    dr[[j]] <- (sp$r - sm$r) / (2 * hj)
    dud[[j]] <- (sp$ud - sm$ud) / (2 * hj)
  }
  uf <- u[base$free]
  rhs_j <- function(j) dr[[j]] - as.numeric(dK[[j]] %*% uf)
  out <- list(solution = base$solution, theta = theta, method = method)
  if (method == "forward") {
    S <- matrix(0, length(u), p)
    for (j in seq_len(p)) {
      S[base$free, j] <- base$solver$solve(rhs_j(j))$x
      if (length(base$dofs)) S[base$dofs, j] <- dud[[j]]
    }
    out$sensitivities <- S
    if (!is.null(objective)) {
      ob <- objective(u)
      out$value <- ob$value
      out$gradient <- as.numeric(crossprod(S, ob$gradient))
    }
  } else {
    ob <- objective(u)
    KT <- if (isTRUE(base$solver$info$symmetric)) base$solver else linear_solver(Matrix::t(base$K))
    lam <- KT$solve(ob$gradient[base$free])$x
    gr <- vapply(seq_len(p), function(j) {
      sum(lam * rhs_j(j)) + if (length(base$dofs)) sum(ob$gradient[base$dofs] * dud[[j]]) else 0
    }, numeric(1))
    out$value <- ob$value
    out$gradient <- gr
    out$adjoint <- lam
  }
  structure(out, class = "pde_sensitivity")
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
  f <- disc$free; d <- disc$dirichlet$dofs
  K <- disc$A[f, f, drop = FALSE]
  r <- disc$b[f]
  if (length(d)) r <- r - as.numeric(disc$A[f, d, drop = FALSE] %*% disc$dirichlet$values)
  out <- list(K = K, r = r, ud = disc$dirichlet$values, free = f, dofs = d)
  if (solve) {
    slv <- linear_solver(K)
    u <- numeric(disc$space$ndof)
    u[d] <- disc$dirichlet$values
    u[f] <- slv$solve(r)$x
    out$u <- u
    out$solver <- slv
    out$solution <- structure(list(u = fe_function(disc$space, u), space = disc$space,
                                   problem = problem), class = "pde_solution")
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
#' @export
pde_observations <- function(x, y, values, sd = 1) {
  n <- length(x)
  if (length(y) != n || length(values) != n) stopf("`x`, `y` and `values` must have equal length.")
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
