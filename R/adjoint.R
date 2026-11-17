#' Adjoint gradients for nonlinear and time-dependent problems
#'
#' `nonlinear_sensitivity()` returns the gradient of `J(u(theta))` for a
#' [nonlinear_problem()] family. At the converged state the discrete
#' residual satisfies `F(u, theta) = 0`; one solve with the transposed Newton
#' Jacobian, `J_F^T lambda = dJ/du`, gives
#' `dJ/dtheta = dJ/du_d du_d/dtheta - lambda' dF/dtheta`, where the
#' parameter derivative of the residual is a central difference at the fixed
#' state (no extra nonlinear solves).
#'
#' `transient_sensitivity()` returns the gradient of an objective of the
#' final state of a [parabolic_problem()] family solved by backward Euler or
#' Crank-Nicolson (without startup steps). The forward run stores the
#' states; a backward sweep with the transposed step matrix computes the
#' discrete adjoint, and the gradient collects the derivatives of the
#' operator, the mass matrix (which depends on the parameters under SUPG
#' stabilization), the loads and the initial condition. It is the exact gradient
#' of the discrete time-stepping scheme (up to the difference steps of the
#' operator derivatives), at the cost of one extra solve per step,
#' independent of the number of parameters. Coefficients must not depend on
#' time and the Dirichlet data must not depend on `theta`. With
#' `checkpoint`, only every `k`-th state is stored and the others are
#' recomputed during the backward sweep.
#'
#' Check either gradient with Taylor remainders, as in [taylor_test()].
#'
#' @param problem_fun `function(theta)` returning a [nonlinear_problem()] or
#'   [parabolic_problem()].
#' @param theta Parameter vector.
#' @param objective `function(u)` of the full degree-of-freedom vector (the
#'   final state for transient problems) returning
#'   `list(value = , gradient = )`, for example from [misfit_objective()].
#' @param h Relative difference step for parameter derivatives.
#' @param dt,method Time step and scheme for transient problems.
#' @param ... Passed to [solve_pde()].
#' @return A list with `value`, `gradient`, the `solution` and the adjoint.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 6)
#' pf <- function(th) nonlinear_problem(mesh, diffusion = function(u) th[1] + u^2, source = th[2],
#'                                      bcs = dirichlet_bc(0))
#' ob <- misfit_objective(fe_space(mesh), pde_observations(0.5, 0.5, 0.05))
#' nonlinear_sensitivity(pf, c(1, 1), ob)$gradient
#' @export
nonlinear_sensitivity <- function(problem_fun, theta, objective, h = 1e-6, ...) {
  theta <- as.numeric(theta); check_tol(h, "h")
  sol <- solve_pde(problem_fun(theta), ...)
  if (!isTRUE(sol$converged)) warnf("Newton did not converge; the gradient is not reliable.")
  ns <- nonlinear_system(problem_fun(theta))
  u <- sol$u$values; f <- ns$free; d <- ns$dd$dofs
  ob <- objective(u)
  J <- ns$residual_jac(u)$J
  as <- solve_linear_system(Matrix::t(J[f, f, drop = FALSE]), ob$gradient[f], method = "lu")
  warn_solver(as$info, "adjoint solve")
  lam <- as$x
  gr <- vapply(seq_along(theta), function(j) {
    hj <- h * max(1, abs(theta[j]))
    res_at <- function(th) {
      nj <- nonlinear_system(problem_fun(th))
      if (!identical(nj$dd$dofs, d)) stopf("The Dirichlet degrees of freedom must not depend on theta.")
      uj <- u; uj[d] <- nj$dd$values
      list(r = nj$residual_jac(uj, jacobian = FALSE)$r, ud = nj$dd$values)
    }
    tp <- theta; tp[j] <- tp[j] + hj; tm <- theta; tm[j] <- tm[j] - hj
    rp <- res_at(tp); rm <- res_at(tm)
    -sum(lam * (rp$r[f] - rm$r[f])) / (2 * hj) +
      if (length(d)) sum(ob$gradient[d] * (rp$ud - rm$ud)) / (2 * hj) else 0
  }, numeric(1))
  list(value = ob$value, gradient = gr, solution = sol, adjoint = lam)
}

#' @rdname nonlinear_sensitivity
#' @param checkpoint For transient problems, `NULL` to store every state of
#'   the forward run, or a step interval `k`: only every `k`-th state is
#'   kept and the states in between are recomputed segment by segment during
#'   the backward sweep (one extra forward solve per step). With
#'   `k = ceiling(sqrt(N))` for `N` steps, memory falls from `N` to about
#'   `2 sqrt(N)` state vectors. The gradient is the same up to rounding.
#' @export
transient_sensitivity <- function(problem_fun, theta, objective, dt,
                                  method = c("backward-euler", "crank-nicolson"), h = 1e-6, checkpoint = NULL, ...) {
  method <- match.arg(method)
  theta <- as.numeric(theta); check_tol(h, "h")
  if (!is.null(checkpoint)) checkpoint <- check_count(checkpoint, "checkpoint", min = 1L)
  th_s <- if (method == "crank-nicolson") 0.5 else 1
  prob <- problem_fun(theta)
  if (!inherits(prob, "parabolic_problem")) stopf("`problem_fun` must return a parabolic_problem().")
  ops <- parabolic_operators(prob)
  if (ops$coef_time_dependent) stopf("Coefficients must not depend on time.")
  V <- ops$V
  t0 <- prob$time[1]; t1 <- prob$time[2]
  N <- max(1L, as.integer(ceiling((t1 - t0) / dt - 1e-10))); dt <- (t1 - t0) / N
  grid <- t0 + (0:N) * dt
  # Forward run: all states, or only every `checkpoint`-th one.
  if (is.null(checkpoint)) {
    sol <- solve_pde(prob, dt = dt, method = method, startup = 0L, ...)
    U <- sol$values
    getU <- function(n) U[, n]
  } else {
    keep <- seq(1L, N + 1L, by = checkpoint)
    ck <- matrix(0, V$ndof, length(keep)); ck[, 1] <- initial_values(prob$initial, V)
    d0 <- dirichlet_data(V, prob$bcs, t0); if (length(d0$dofs)) ck[d0$dofs, 1] <- d0$values
    step <- 0L
    sol <- solve_pde(prob, dt = dt, method = method, startup = 0L, output_times = t1, ...,
                     callback = function(t, u) {
                       step <<- step + 1L
                       k <- match(step + 1L, keep); if (!is.na(k)) ck[, k] <<- u
                     })
    seg <- NULL; seg_lo <- 0L
    getU <- function(n) {
      if (is.null(seg) || n < seg_lo || n > seg_lo + ncol(seg) - 1L) {
        c0 <- keep[findInterval(n, keep)]
        c1 <- min(c0 + checkpoint, N + 1L)
        if (c1 == c0) { seg <<- ck[, match(c0, keep), drop = FALSE]; seg_lo <<- c0; return(seg[, 1]) }
        pr <- prob; pr$time <- c(t0, grid[c1])
        rs <- structure(list(state = list(kind = "parabolic", t = grid[c0], u = ck[, match(c0, keep)],
                                          step = c0 - 1L, method = method, version = state_version, space = space_id(V))),
                       class = "pde_transient")
        sg <- solve_pde(pr, dt = dt, method = method, startup = 0L, restart = rs, ...)
        seg <<- sg$values; seg_lo <<- c0
      }
      seg[, n - seg_lo + 1L]
    }
  }
  uN <- sol$state$u
  d <- dirichlet_data(V, prob$bcs, t0)$dofs
  f <- setdiff(seq_len(V$ndof), d)
  L <- ops$M + (th_s * dt) * ops$A; R <- ops$M - ((1 - th_s) * dt) * ops$A
  Lt <- linear_solver(Matrix::t(L[f, f, drop = FALSE]))
  Rff_t <- Matrix::t(R[f, f, drop = FALSE])
  ob <- objective(uN)
  # Parameter derivatives of the operator, loads and initial state.
  der <- lapply(seq_along(theta), function(j) {
    hj <- h * max(1, abs(theta[j]))
    tp <- theta; tp[j] <- tp[j] + hj; tm <- theta; tm[j] <- tm[j] - hj
    pp <- problem_fun(tp); pm <- problem_fun(tm)
    # The adjoint omits derivatives of prescribed values; check that there
    # are none (at the first, middle and last time).
    for (t in grid[unique(c(1L, N %/% 2L + 1L, N + 1L))]) {
      dp <- dirichlet_data(V, pp$bcs, t); dm <- dirichlet_data(V, pm$bcs, t)
      if (!identical(dp$dofs, d) || !identical(dm$dofs, d) ||
          max(abs(dp$values - dm$values), 0) > 1e-12 * max(1, abs(dp$values))) {
        stopf("The Dirichlet data must not depend on theta (parameter %d).", j)
      }
    }
    op <- parabolic_operators(pp); om <- parabolic_operators(pm)
    # SUPG makes the mass matrix depend on the velocity and diffusion, so
    # its derivative enters the step residual too.
    list(dA = (op$A - om$A) / (2 * hj), dM = (op$M - om$M) / (2 * hj),
         du0 = c_zero((initial_values(pp$initial, V) - initial_values(pm$initial, V)) / (2 * hj), d),
         dF = function(t) (op$load(t) - om$load(t)) / (2 * hj))
  })
  # Backward sweep: L' lambda_N = g, L' lambda_n = R' lambda_{n+1}; each
  # step's contribution to every gradient entry is added as lambda_n is
  # formed, so neither the adjoints nor (with checkpoints) the states are
  # stored.
  gr <- numeric(length(theta))
  dF_next <- lapply(der, function(dj) dj$dF(grid[N + 1L]))
  u_next <- uN
  lam <- NULL; full_lam <- if (is.null(checkpoint)) matrix(0, length(f), N) else NULL
  for (n in N:1L) {
    s <- Lt$solve(if (n == N) ob$gradient[f] else as.numeric(Rff_t %*% lam))
    warn_solver(s$info, "adjoint solve")
    lam <- s$x
    if (!is.null(full_lam)) full_lam[, n] <- lam
    u_n <- getU(n)
    for (j in seq_along(der)) {
      dF_n <- der[[j]]$dF(grid[n])
      rhs <- -((1 - th_s) * dt) * as.numeric(der[[j]]$dA %*% u_n) - (th_s * dt) * as.numeric(der[[j]]$dA %*% u_next) +
        dt * (th_s * dF_next[[j]] + (1 - th_s) * dF_n) - as.numeric(der[[j]]$dM %*% (u_next - u_n))
      gr[j] <- gr[j] + sum(lam * rhs[f])
      dF_next[[j]] <- dF_n
    }
    u_next <- u_n
  }
  for (j in seq_along(der)) gr[j] <- gr[j] + sum(lam * as.numeric(R %*% der[[j]]$du0)[f])
  # (The step residual is L u_n - R u_{n-1} - dt F with L = M + theta dt A
  # and R = M - (1 - theta) dt A, so dM enters as dM (u_n - u_{n-1}).)
  list(value = ob$value, gradient = gr, solution = sol, adjoint = full_lam)
}

# Zero the prescribed entries (their values do not depend on theta).
c_zero <- function(v, d) { v[d] <- 0; v }
