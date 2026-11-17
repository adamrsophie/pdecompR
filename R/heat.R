#' Linear parabolic problems
#'
#' `parabolic_problem()` defines
#'
#' `du/dt - div(k grad u) + b . grad u + c u = f(x, y, t)` for `t` in
#' `(time[1], time[2])`, with `u(x, y, time[1]) = initial(x, y)` and
#' boundary conditions whose data may depend on time (`function(x, y, t)`).
#' Scalar diffusion, advection and reaction coefficients may also depend on
#' time (`function(x, y, t)`); the operator is then reassembled and
#' refactorized at every step, with `A(t_{n+1})` on the implicit side and
#' `A(t_n)` on the explicit side of the theta scheme (second order for
#' Crank-Nicolson). Otherwise the matrices are assembled and factorized once.
#' The Rannacher startup steps also apply with time-dependent coefficients.
#' SUPG requires time-independent coefficients.
#'
#' `heat_problem()` is the heat equation `du/dt - kappa Laplace(u) = f`.
#'
#' The finite-element semidiscretization is `M du/dt + A u = F(t)` with mass
#' matrix `M` and operator matrix `A` (including Robin terms). Time stepping
#' uses the theta method; see [solve_pde()] arguments `dt` and `method`.
#'
#' @inheritParams elliptic_problem
#' @param initial Initial condition: number, `function(x, y)` or
#'   [fe_function()].
#' @param time Start and end time.
#' @param diffusivity Diffusivity `kappa` of the heat equation.
#' @return An object of class `parabolic_problem`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' prob <- heat_problem(mesh, diffusivity = 1,
#'                      initial = function(x, y) sin(pi * x) * sin(pi * y),
#'                      bcs = dirichlet_bc(0), time = c(0, 0.1))
#' sol <- solve_pde(prob, dt = 0.01)
#' sol
#' @export
parabolic_problem <- function(mesh, diffusion = 1, reaction = 0, advection = NULL, source = 0,
                              initial = 0, bcs = list(), time = c(0, 1), degree = 1,
                              stabilization = c("none", "supg"), quad_degree = NULL) {
  ep <- elliptic_problem(mesh, diffusion = diffusion, reaction = reaction, advection = advection,
                         source = 0, bcs = bcs, degree = degree,
                         stabilization = match.arg(stabilization), quad_degree = quad_degree)
  if (!is.numeric(time) || length(time) != 2L || !all(is.finite(time)) || time[2] <= time[1]) {
    stopf("`time` must be two increasing finite numbers.")
  }
  if (!(is.function(source) || (is.numeric(source) && length(source) == 1L))) {
    stopf("`source` must be a number or a function of (x, y) or (x, y, t).")
  }
  ep$source <- source
  ep$initial <- initial
  ep$time <- as.numeric(time)
  class(ep) <- c("parabolic_problem", "pde_problem")
  ep
}

#' @rdname parabolic_problem
#' @export
heat_problem <- function(mesh, diffusivity = 1, source = 0, initial = 0, bcs = list(),
                         time = c(0, 1), degree = 1) {
  parabolic_problem(mesh, diffusion = diffusivity, source = source, initial = initial,
                    bcs = bcs, time = time, degree = degree)
}

# Time-independent matrices and a load function F(t) for a parabolic problem.
parabolic_operators <- function(problem) {
  t0 <- proc.time()[["elapsed"]]
  V <- fe_space(problem$mesh, problem$degree)
  geo <- element_geometry(V, problem$quad_degree)
  coef_td <- is_time_function(problem$diffusion) || is_time_function(problem$reaction) ||
    is_time_function(problem$advection)
  if (coef_td && problem$stabilization == "supg") {
    stopf("SUPG with time-dependent coefficients is not supported.")
  }
  t_init <- problem$time[1]
  D <- diffusion_qp(problem$diffusion, geo, t_init)
  if (!is.list(D) && min(D) <= 0) stopf("The diffusion coefficient must be positive.")
  cq <- coef_qp(problem$reaction, geo, t_init, name = "reaction")
  if (all(cq == 0)) cq <- NULL
  bq <- vec_qp(problem$advection, geo, t_init)
  tq <- NULL
  advt <- NULL; madvt <- NULL
  if (problem$stabilization == "supg") {
    tq <- matrix(supg_tau(problem$mesh, problem$diffusion, problem$advection), geo$m, geo$nq)
    D <- supg_diffusion(D, tq, bq)
    if (!is.null(cq)) advt <- list(x = tq * cq * bq$x, y = tq * cq * bq$y)
    madvt <- list(x = tq * bq$x, y = tq * bq$y)
  }
  M <- global_matrix(V, local_bilinear(geo, mass = matrix(1, geo$m, geo$nq), advt = madvt))
  A <- global_matrix(V, local_bilinear(geo, mass = cq, adv = bq, advt = advt, diff = D))
  rb <- Filter(function(b) b$type == "robin", problem$bcs)
  if (any(vapply(rb, function(b) is_time_function(b$alpha), logical(1)))) {
    stopf("Time-dependent Robin coefficients `alpha` are not supported; the Robin value may depend on t.")
  }
  robin <- natural_bc_terms(V, rb, problem$time[1])$A
  if (!is.null(robin)) A <- A + robin
  # Operator at time t, reassembled when a coefficient depends on time.
  A_at <- function(t) {
    if (!coef_td) return(A)
    Dt <- diffusion_qp(problem$diffusion, geo, t)
    if (!is.list(Dt) && min(Dt) <= 0) stopf("The diffusion coefficient must be positive (t = %s).", fmt_num(t))
    ct <- coef_qp(problem$reaction, geo, t, name = "reaction")
    if (all(ct == 0)) ct <- NULL
    At <- global_matrix(V, local_bilinear(geo, mass = ct, adv = vec_qp(problem$advection, geo, t), diff = Dt))
    if (!is.null(robin)) At <- At + robin
    At
  }
  load <- function(t) {
    fq <- coef_qp(problem$source, geo, t, "source")
    g <- if (!is.null(tq)) list(x = tq * fq * bq$x, y = tq * fq * bq$y) else NULL
    global_vector(V, local_linear(geo, f = fq, g = g)) +
      natural_bc_terms(V, problem$bcs, t, problem$quad_degree, matrix = FALSE)$b
  }
  time_dep <- is_time_function(problem$source) ||
    any(vapply(problem$bcs, function(b) is_time_function(b$value), logical(1)))
  list(V = V, M = M, A = A, A_at = A_at, coef_time_dependent = coef_td, load = load,
       time_dependent = time_dep,
       assembly_time = proc.time()[["elapsed"]] - t0)
}

#' @rdname solve_pde
#' @param dt Time step. The number of steps is `ceiling((T - t0) / dt)` and
#'   the step is then adjusted to land exactly on the final time.
#' @param method For time-dependent problems: `"crank-nicolson"` (second
#'   order), `"backward-euler"` (first order, L-stable), `"bdf2"` (second
#'   order, L-stable two-step method started with one backward Euler step) or
#'   `"sdirk2"` (Alexander's two-stage, second-order, L-stable singly
#'   diagonally implicit Runge-Kutta method). BDF2 and SDIRK2 require
#'   coefficients that do not depend on time (data may).
#' @param output_times Times at which to store the solution; they must lie on
#'   the time grid. Default: every step. With `adaptive = TRUE` they may be
#'   any times in the interval, and the steps land on them exactly.
#' @param adaptive For `"sdirk2"`: choose the step size from an embedded
#'   first-order error estimate (L2 norm through the mass matrix), starting
#'   from `dt`. Steps whose estimate exceeds `tol` are rejected and retried.
#' @param tol Error tolerance per step for adaptive stepping, relative to
#'   `max(1, ||u||)`.
#' @param solver Linear solver method, see [solve_linear_system()].
#' @param startup Number of initial backward Euler steps taken before
#'   Crank-Nicolson (Rannacher startup). They damp the oscillations that
#'   Crank-Nicolson produces for rough initial data or initial data that do
#'   not match the Dirichlet data; a fixed number of such steps keeps second
#'   order. Ignored for backward Euler.
#' @references Rannacher, R. (1984). Finite element solution of diffusion
#'   problems with irregular data. *Numerische Mathematik*, 43, 309-327.
#' @export
solve_pde.parabolic_problem <- function(problem, dt, method = c("crank-nicolson", "backward-euler", "bdf2", "sdirk2"),
                                        output_times = NULL, solver = "auto", startup = 2L,
                                        adaptive = FALSE, tol = 1e-4, ...) {
  method <- match.arg(method)
  check_number(dt, "dt", positive = TRUE)
  if (!is.logical(adaptive) || length(adaptive) != 1L || is.na(adaptive)) stopf("`adaptive` must be TRUE or FALSE.")
  if (adaptive && method != "sdirk2") stopf("Adaptive time stepping is available for method = \"sdirk2\".")
  if (method %in% c("bdf2", "sdirk2")) {
    return(solve_parabolic_multistage(problem, dt, method, output_times, solver, adaptive, tol, ...))
  }
  startup <- check_count(startup, "startup", min = 0L)
  theta <- if (method == "crank-nicolson") 0.5 else 1
  if (theta == 1) startup <- 0L
  ops <- parabolic_operators(problem)
  V <- ops$V
  t0 <- problem$time[1]; t1 <- problem$time[2]
  N <- max(1L, as.integer(ceiling((t1 - t0) / dt - 1e-10)))
  dt <- (t1 - t0) / N
  startup <- min(startup, N)
  grid <- t0 + (0:N) * dt
  out_idx <- output_index(output_times, grid, dt)
  dd0 <- dirichlet_data(V, problem$bcs, t0)
  bdofs <- dd0$dofs
  free <- setdiff(seq_len(V$ndof), bdofs)
  u <- initial_values(problem$initial, V)
  if (length(bdofs)) u[bdofs] <- dd0$values
  L <- ops$M + (theta * dt) * ops$A
  R <- ops$M - ((1 - theta) * dt) * ops$A
  ts <- proc.time()[["elapsed"]]
  slv <- linear_solver(L[free, free, drop = FALSE], solver, ...)
  Lfd <- L[free, bdofs, drop = FALSE]
  A_old <- ops$A
  if (startup > 0L && !ops$coef_time_dependent) {
    Lbe <- ops$M + dt * ops$A
    slv_be <- linear_solver(Lbe[free, free, drop = FALSE], solver, ...)
    Lfd_be <- Lbe[free, bdofs, drop = FALSE]
  }
  U <- matrix(NA_real_, V$ndof, length(out_idx))
  k_out <- 1L
  if (out_idx[1] == 1L) { U[, 1] <- u; k_out <- 2L }
  F_old <- ops$load(t0)
  max_rel <- 0
  bad <- NULL
  for (n in seq_len(N)) {
    tn <- grid[n + 1L]
    F_new <- if (ops$time_dependent) ops$load(tn) else F_old
    be <- n <= startup
    if (ops$coef_time_dependent) {
      # Theta scheme with A(t_{n+1}) on the left and A(t_n) on the right.
      A_new <- ops$A_at(tn)
      th <- if (be) 1 else theta
      Ln <- ops$M + (th * dt) * A_new
      slv <- linear_solver(Ln[free, free, drop = FALSE], solver, ...)
      Lfd <- Ln[free, bdofs, drop = FALSE]
      R <- ops$M - ((1 - th) * dt) * A_old
      A_old <- A_new
      be <- FALSE
    }
    rhs <- if (be) as.numeric(ops$M %*% u) + dt * F_new
           else if (ops$coef_time_dependent) as.numeric(R %*% u) + dt * (th * F_new + (1 - th) * F_old)
           else as.numeric(R %*% u) + dt * (theta * F_new + (1 - theta) * F_old)
    rhs <- rhs[free]
    ud <- if (length(bdofs)) dirichlet_data(V, problem$bcs, tn)$values else numeric()
    if (length(bdofs)) rhs <- rhs - as.numeric((if (be) Lfd_be else Lfd) %*% ud)
    s <- (if (be) slv_be else slv)$solve(rhs, x0 = u[free])  # warm start for iterative methods
    u[free] <- s$x
    if (length(bdofs)) u[bdofs] <- ud
    max_rel <- max(max_rel, s$info$relative_residual)
    if (!s$info$status %in% c("solved", "converged")) bad <- s$info
    F_old <- F_new
    if (k_out <= length(out_idx) && out_idx[k_out] == n + 1L) {
      U[, k_out] <- u; k_out <- k_out + 1L
    }
  }
  if (!is.null(bad)) warn_solver(bad, "linear solve of at least one time step")
  diag <- structure(list(ndof = V$ndof, free = length(free), steps = N, dt = dt,
                         method = method, theta = theta, startup = startup, solver = slv$info,
                         integral = as.numeric(crossprod(assemble_vector(V, 1), U)),
                         max_relative_residual = max_rel,
                         assembly_time = ops$assembly_time,
                         solve_time = proc.time()[["elapsed"]] - ts),
                    class = "pde_diagnostics")
  structure(list(times = grid[out_idx], values = U, space = V, problem = problem,
                 diagnostics = diag), class = "pde_transient")
}

output_index <- function(output_times, grid, dt) {
  if (is.null(output_times)) return(seq_along(grid))
  if (!is.numeric(output_times) || !length(output_times) || any(!is.finite(output_times))) {
    stopf("`output_times` must be a nonempty vector of finite times.")
  }
  k <- round((output_times - grid[1]) / dt)
  if (any(k < 0 | k > length(grid) - 1L) ||
      any(abs(grid[1] + k * dt - output_times) > 1e-8 * max(1, abs(output_times)))) {
    stopf("`output_times` must lie on the time grid t0 + k * dt (dt = %s).", fmt_num(dt, 8))
  }
  sort(unique(k + 1L))
}

initial_values <- function(initial, V) {
  if (inherits(initial, "fe_function")) {
    if (same_space(initial$space, V)) {
      return(initial$values)
    }
    return(transfer_function(initial, V$mesh, V$degree)$values)
  }
  if (is.numeric(initial) && length(initial) == V$ndof && V$ndof > 1L) return(as.numeric(initial))
  eval_spatial(initial, V$dof_coords[, 1], V$dof_coords[, 2], name = "initial")
}

#' Extract a time snapshot
#'
#' @param solution A `pde_transient` solution.
#' @param time A stored output time (the nearest stored time is used).
#' @param species For reaction-diffusion systems, the species name or index.
#' @return An [fe_function()].
#' @examples
#' prob <- heat_problem(pde_mesh_rectangle(nx = 4), initial = function(x, y) sin(pi * x) * sin(pi * y),
#'                      bcs = dirichlet_bc(0), time = c(0, 0.1))
#' sol <- solve_pde(prob, dt = 0.05)
#' u <- solution_at(sol, 0.1)
#' evaluate(u, 0.5, 0.5)
#' @export
solution_at <- function(solution, time, species = 1) {
  stopifnot(inherits(solution, "pde_transient"))
  k <- which.min(abs(solution$times - time))
  if (abs(solution$times[k] - time) > 1e-8 * max(1, abs(time))) {
    warnf("Time %s is not stored; using the nearest stored time %s.",
          fmt_num(time), fmt_num(solution$times[k]))
  }
  vals <- if (is.list(solution$values)) solution$values[[species]][, k] else solution$values[, k]
  nm <- if (is.list(solution$values)) {
    if (is.character(species)) species else names(solution$values)[species]
  } else "u"
  fe_function(solution$space, vals, name = nm)
}

#' @export
print.pde_transient <- function(x, ...) {
  d <- x$diagnostics
  cat("<pde_transient> P", x$space$degree, ", ", d$ndof, " dofs, ", d$steps,
      " steps of ", d$method, " (dt = ", fmt_num(d$dt), ")\n", sep = "")
  cat("  stored times: ", length(x$times), " in [", fmt_num(min(x$times)), ", ",
      fmt_num(max(x$times)), "]\n", sep = "")
  invisible(x)
}

#' @export
summary.pde_transient <- function(object, ...) {
  d <- object$diagnostics
  cat("Time-dependent finite-element solution\n")
  cat("  method:", d$method, " dt:", fmt_num(d$dt), " steps:", d$steps, "\n")
  cat("  degrees of freedom:", d$ndof, " (", d$free, " free)\n")
  cat("  linear solver:", d$solver$method, " max relative residual:",
      fmt_num(d$max_relative_residual, 3), "\n")
  if (!is.null(d$newton_iterations)) {
    cat("  Newton iterations per step: max", max(d$newton_iterations), " mean",
        fmt_num(mean(d$newton_iterations)), "\n")
  }
  if (!is.null(d$integral)) {
    I <- as.matrix(d$integral)
    if (nrow(I) == 1L && ncol(I) > 1L && is.null(colnames(I))) I <- t(I)
    cat("  integral of the solution, first and last stored time:",
        paste(fmt_num(I[1, ]), "to", fmt_num(I[nrow(I), ]), collapse = "; "), "\n")
  }
  cat("  assembly", fmt_num(d$assembly_time, 3), "s, time stepping", fmt_num(d$solve_time, 3), "s\n")
  invisible(object)
}

# BDF2 and SDIRK2 (fixed or adaptive step) for M u' + A u = F(t) with
# Dirichlet data g(t). Stage equations are written for the stage values U_i
# with their prescribed boundary values, as in the theta scheme, so no time
# derivative of the Dirichlet data is needed.
solve_parabolic_multistage <- function(problem, dt, method, output_times, solver, adaptive, tol, ...) {
  check_tol(tol, "tol")
  ops <- parabolic_operators(problem)
  if (ops$coef_time_dependent) {
    stopf("Method \"%s\" requires coefficients that do not depend on time; use crank-nicolson.", method)
  }
  V <- ops$V; M <- ops$M; A <- ops$A
  t0 <- problem$time[1]; t1 <- problem$time[2]
  dd0 <- dirichlet_data(V, problem$bcs, t0)
  bdofs <- dd0$dofs
  free <- setdiff(seq_len(V$ndof), bdofs)
  u <- initial_values(problem$initial, V)
  if (length(bdofs)) u[bdofs] <- dd0$values
  load <- function(t) if (ops$time_dependent) ops$load(t) else F0
  F0 <- ops$load(t0)
  gd <- function(t) if (length(bdofs)) dirichlet_data(V, problem$bcs, t)$values else numeric()
  ts <- proc.time()[["elapsed"]]
  max_rel <- 0; bad <- NULL
  note <- function(s) {
    max_rel <<- max(max_rel, s$info$relative_residual)
    if (!solve_ok(s$info)) bad <<- s$info
  }
  # Solve (c M + h A) U = rhs on the free rows with U = g on bdofs.
  cache <- new.env()
  stage_solve <- function(cm, h, rhs, g, x0) {
    key <- sprintf("%.17g/%.17g", cm, h)
    if (!identical(cache$key, key)) {
      L <- cm * M + h * A
      cache$key <- key
      cache$slv <- linear_solver(L[free, free, drop = FALSE], solver, ...)
      cache$Lfd <- L[free, bdofs, drop = FALSE]
      cache$n <- (if (is.null(cache$n)) 0L else cache$n) + 1L
    }
    r <- rhs[free]
    if (length(bdofs)) r <- r - as.numeric(cache$Lfd %*% g)
    s <- cache$slv$solve(r, x0 = x0[free]); note(s)
    U <- numeric(V$ndof); U[free] <- s$x
    if (length(bdofs)) U[bdofs] <- g
    U
  }
  gam <- 1 - 1 / sqrt(2)
  sdirk <- function(u, t, h) {
    # Alexander's SDIRK2: A = [[g, 0], [1 - g, g]], c = (g, 1), stiffly accurate.
    F1 <- load(t + gam * h)
    U1 <- stage_solve(1, gam * h, as.numeric(M %*% u) + gam * h * F1, gd(t + gam * h), u)
    Y1 <- F1 - as.numeric(A %*% U1)
    F2 <- load(t + h)
    U2 <- stage_solve(1, gam * h, as.numeric(M %*% u) + (1 - gam) * h * Y1 + gam * h * F2, gd(t + h), U1)
    # Embedded first-order solution u + h K1 with K1 = (U1 - u) / (g h).
    uh <- u + (U1 - u) / gam
    list(u = U2, err = U2 - uh)
  }
  mnorm <- function(v) sqrt(max(0, sum(v * as.numeric(M %*% v))))
  if (!adaptive) {
    N <- max(1L, as.integer(ceiling((t1 - t0) / dt - 1e-10)))
    dt <- (t1 - t0) / N
    grid <- t0 + (0:N) * dt
    out_idx <- output_index(output_times, grid, dt)
    U <- matrix(NA_real_, V$ndof, length(out_idx))
    k_out <- 1L
    if (out_idx[1] == 1L) { U[, 1] <- u; k_out <- 2L }
    u_old <- NULL
    for (n in seq_len(N)) {
      tn <- grid[n + 1L]
      u_new <- if (method == "sdirk2") sdirk(u, grid[n], dt)$u
               else if (is.null(u_old)) stage_solve(1, dt, as.numeric(M %*% u) + dt * load(tn), gd(tn), u)
               else stage_solve(1.5, dt, as.numeric(M %*% (2 * u - 0.5 * u_old)) + dt * load(tn), gd(tn), u)
      u_old <- u; u <- u_new
      if (k_out <= length(out_idx) && out_idx[k_out] == n + 1L) { U[, k_out] <- u; k_out <- k_out + 1L }
    }
    times <- grid[out_idx]; steps <- N; rejected <- 0L; dts <- rep(dt, N)
  } else {
    outs <- if (is.null(output_times)) NULL else {
      if (!is.numeric(output_times) || !length(output_times) || any(!is.finite(output_times)) ||
          any(output_times < t0 - 1e-12 * max(1, abs(t0)) | output_times > t1 + 1e-12 * max(1, abs(t1)))) {
        stopf("`output_times` must be finite times in [%s, %s].", fmt_num(t0), fmt_num(t1))
      }
      sort(unique(output_times))
    }
    stops <- sort(unique(c(outs, t1)))
    tt <- t0; h <- min(dt, t1 - t0)
    times <- t0; vals <- list(u); dts <- numeric(0); rejected <- 0L
    hmin <- 1e-12 * max(1, abs(t1 - t0))
    while (tt < t1 - hmin) {
      nxt <- stops[stops > tt + hmin][1]
      hs <- min(h, nxt - tt)
      st <- sdirk(u, tt, hs)
      err <- mnorm(st$err) / max(1, mnorm(st$u))
      if (!is.finite(err)) stop_nonfinite("The time step produced non-finite values at t = %s.", fmt_num(tt))
      fac <- min(2, max(0.2, 0.9 * sqrt(tol / max(err, 1e-300))))
      if (err <= tol) {
        tt <- if (abs(tt + hs - nxt) <= hmin) nxt else tt + hs
        u <- st$u; dts <- c(dts, hs)
        if (is.null(outs) || any(abs(outs - tt) <= hmin)) { times <- c(times, tt); vals[[length(vals) + 1L]] <- u }
        h <- hs * fac
      } else {
        rejected <- rejected + 1L
        h <- hs * fac
        if (h < hmin) stopf("The adaptive step size fell below %s at t = %s.", fmt_num(hmin), fmt_num(tt))
      }
    }
    keep <- if (is.null(outs)) rep(TRUE, length(times)) else vapply(times, function(t) any(abs(outs - t) <= hmin), logical(1))
    times <- times[keep]; U <- do.call(cbind, vals[keep]); steps <- length(dts)
  }
  if (!is.null(bad)) warn_solver(bad, "linear solve of at least one time step")
  diag <- structure(list(ndof = V$ndof, free = length(free), steps = steps,
                         dt = if (adaptive) range(dts) else dt, method = method, adaptive = adaptive,
                         rejected_steps = rejected, factorizations = cache$n,
                         integral = as.numeric(crossprod(assemble_vector(V, 1), U)),
                         max_relative_residual = max_rel, assembly_time = ops$assembly_time,
                         solve_time = proc.time()[["elapsed"]] - ts),
                    class = "pde_diagnostics")
  structure(list(times = times, values = U, space = V, problem = problem, diagnostics = diag),
            class = "pde_transient")
}
