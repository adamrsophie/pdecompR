#' Linear parabolic problems
#'
#' `parabolic_problem()` defines
#'
#' `du/dt - div(k grad u) + b . grad u + c u = f(x, y, t)` for `t` in
#' `(time[1], time[2])`, with `u(x, y, time[1]) = initial(x, y)` and
#' boundary conditions whose data may depend on time (`function(x, y, t)`).
#' The diffusion, advection and reaction coefficients are time independent,
#' so the system matrices are assembled once.
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
  D <- diffusion_qp(problem$diffusion, geo)
  if (!is.list(D) && min(D) <= 0) stopf("The diffusion coefficient must be positive.")
  cq <- coef_qp(problem$reaction, geo, name = "reaction")
  if (all(cq == 0)) cq <- NULL
  bq <- vec_qp(problem$advection, geo)
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
  load <- function(t) {
    fq <- coef_qp(problem$source, geo, t, "source")
    g <- if (!is.null(tq)) list(x = tq * fq * bq$x, y = tq * fq * bq$y) else NULL
    global_vector(V, local_linear(geo, f = fq, g = g)) +
      natural_bc_terms(V, problem$bcs, t, problem$quad_degree)$b
  }
  time_dep <- is_time_function(problem$source) ||
    any(vapply(problem$bcs, function(b) is_time_function(b$value), logical(1)))
  list(V = V, M = M, A = A, load = load, time_dependent = time_dep,
       assembly_time = proc.time()[["elapsed"]] - t0)
}

#' @rdname solve_pde
#' @param dt Time step. The number of steps is `ceiling((T - t0) / dt)` and
#'   the step is then adjusted to land exactly on the final time.
#' @param method For time-dependent problems: `"crank-nicolson"` (second
#'   order) or `"backward-euler"` (first order, L-stable).
#' @param output_times Times at which to store the solution; they must lie on
#'   the time grid. Default: every step.
#' @param solver Linear solver method, see [solve_linear_system()].
#' @param startup Number of initial backward Euler steps taken before
#'   Crank-Nicolson (Rannacher startup). They damp the oscillations that
#'   Crank-Nicolson produces for rough initial data or initial data that do
#'   not match the Dirichlet data; a fixed number of such steps keeps second
#'   order. Ignored for backward Euler.
#' @references Rannacher, R. (1984). Finite element solution of diffusion
#'   problems with irregular data. *Numerische Mathematik*, 43, 309-327.
#' @export
solve_pde.parabolic_problem <- function(problem, dt, method = c("crank-nicolson", "backward-euler"),
                                        output_times = NULL, solver = "auto", startup = 2L, ...) {
  method <- match.arg(method)
  check_number(dt, "dt", positive = TRUE)
  startup <- check_count(startup, "startup", min = 0L)
  theta <- if (method == "crank-nicolson") 0.5 else 1
  if (theta == 1) startup <- 0L
  ops <- parabolic_operators(problem)
  V <- ops$V
  t0 <- problem$time[1]; t1 <- problem$time[2]
  N <- max(1L, as.integer(ceiling((t1 - t0) / dt - 1e-10)))
  dt <- (t1 - t0) / N
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
  if (startup > 0L) {
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
    rhs <- if (be) as.numeric(ops$M %*% u) + dt * F_new
           else as.numeric(R %*% u) + dt * (theta * F_new + (1 - theta) * F_old)
    rhs <- rhs[free]
    ud <- if (length(bdofs)) dirichlet_data(V, problem$bcs, tn)$values else numeric()
    if (length(bdofs)) rhs <- rhs - as.numeric((if (be) Lfd_be else Lfd) %*% ud)
    s <- (if (be) slv_be else slv)$solve(rhs)
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
                         max_relative_residual = max_rel,
                         assembly_time = ops$assembly_time,
                         solve_time = proc.time()[["elapsed"]] - ts),
                    class = "pde_diagnostics")
  structure(list(times = grid[out_idx], values = U, space = V, problem = problem,
                 diagnostics = diag), class = "pde_transient")
}

output_index <- function(output_times, grid, dt) {
  if (is.null(output_times)) return(seq_along(grid))
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
  cat("  assembly", fmt_num(d$assembly_time, 3), "s, time stepping", fmt_num(d$solve_time, 3), "s\n")
  invisible(object)
}
