#' Wave equation and elastodynamics
#'
#' `wave_problem()` defines the damped wave equation
#'
#' `d2u/dt2 + gamma du/dt - div(c^2 grad u) = f(x, y, t)`
#'
#' with initial displacement and velocity and Dirichlet, Neumann or Robin
#' conditions (data may depend on time). `elastodynamics_problem()` defines
#' `rho d2u/dt2 + gamma rho du/dt - div sigma(u) = f` for the displacement
#' of an [elasticity_problem()], with the static loads of that problem
#' multiplied by `load_factor(t)`.
#'
#' Both semidiscretizations `M a + C v + K u = F(t)` are integrated by the
#' Newmark average acceleration method (`beta = 1/4`, `gamma = 1/2`): second
#' order, unconditionally stable and, without damping and forcing, exactly
#' conserving the discrete energy `1/2 v'Mv + 1/2 u'Ku`, which is reported as
#' `diagnostics$energy`. The effective matrix is factorized once.
#'
#' @param mesh A [pde_mesh()].
#' @param speed Wave speed `c` (number or `function(x, y)`).
#' @param damping Damping coefficient `gamma >= 0`.
#' @param source Source `f`, a number or `function(x, y)` or `function(x, y, t)`.
#' @param initial,velocity Initial displacement and velocity (numbers,
#'   functions of `(x, y)` or finite-element functions).
#' @param bcs Boundary conditions as for [elliptic_problem()].
#' @param time Start and end time.
#' @param degree Polynomial degree.
#' @return An object of class `wave_problem` or `elastodynamics_problem`; solve
#'   it with [solve_pde()] and a time step `dt`.
#' @references Newmark, N. M. (1959). A method of computation for structural
#'   dynamics. *Journal of the Engineering Mechanics Division*, 85, 67-94.
#' @examples
#' u0 <- function(x, y) sin(pi * x) * sin(pi * y)
#' p <- wave_problem(pde_mesh_rectangle(nx = 8), initial = u0, bcs = dirichlet_bc(0),
#'                   time = c(0, 1))
#' s <- solve_pde(p, dt = 0.02)
#' range(s$diagnostics$energy)
#' @export
wave_problem <- function(mesh, speed = 1, damping = 0, source = 0, initial = 0, velocity = 0,
                         bcs = list(), time = c(0, 1), degree = 1) {
  check_coefficient(speed, "speed")
  if (!is.numeric(damping) || length(damping) != 1L || !is.finite(damping) || damping < 0) {
    stopf("`damping` must be a nonnegative number.")
  }
  if (is_time_function(speed)) stopf("`speed` must not depend on time.")
  c2 <- if (is.function(speed)) function(x, y) speed(x, y)^2 else speed^2
  p <- parabolic_problem(mesh, diffusion = c2, source = source, initial = initial, bcs = bcs,
                         time = time, degree = degree)
  p$speed <- speed; p$damping <- damping; p$velocity <- velocity
  class(p) <- c("wave_problem", "pde_problem")
  p
}

#' @rdname wave_problem
#' @param elasticity An [elasticity_problem()] giving the material, loads
#'   and boundary conditions (displacement data are held fixed in time).
#' @param density Mass density `rho > 0` (a number).
#' @param load_factor `function(t)` scaling the body force and tractions.
#' @export
elastodynamics_problem <- function(elasticity, density = 1, damping = 0, initial = c(0, 0),
                                   velocity = c(0, 0), time = c(0, 1), load_factor = function(t) 1) {
  stopifnot(inherits(elasticity, "elasticity_problem"))
  check_number(density, "density", positive = TRUE)
  if (!is.numeric(damping) || length(damping) != 1L || !is.finite(damping) || damping < 0) {
    stopf("`damping` must be a nonnegative number.")
  }
  if (!is.function(load_factor)) stopf("`load_factor` must be a function of t.")
  if (!is.numeric(time) || length(time) != 2L || !all(is.finite(time)) || time[2] <= time[1]) {
    stopf("`time` must be two increasing finite numbers.")
  }
  structure(list(elasticity = elasticity, mesh = elasticity$mesh, degree = elasticity$degree,
                 density = density, damping = damping, initial = initial, velocity = velocity,
                 time = as.numeric(time), load_factor = load_factor),
            class = c("elastodynamics_problem", "pde_problem"))
}

#' @rdname solve_pde
#' @export
solve_pde.wave_problem <- function(problem, dt, output_times = NULL, solver = "auto", ...) {
  ops <- parabolic_operators(problem)
  V <- ops$V
  dd <- function(t) dirichlet_data(V, problem$bcs, t)
  bdofs <- dd(problem$time[1])$dofs
  F0 <- ops$load(problem$time[1])
  load <- if (ops$time_dependent) ops$load else function(t) F0
  out <- newmark(ops$M, problem$damping * ops$M, ops$A, load, bdofs,
                 function(t) dd(t)$values, initial_values(problem$initial, V),
                 initial_values(problem$velocity, V), problem$time, dt, output_times, solver, ...)
  out$diagnostics$assembly_time <- ops$assembly_time
  structure(list(times = out$times, values = out$U, velocity = out$Vel, space = V, problem = problem,
                 diagnostics = out$diagnostics), class = "pde_transient")
}

#' @rdname solve_pde
#' @export
solve_pde.elastodynamics_problem <- function(problem, dt, output_times = NULL, solver = "auto", ...) {
  disc <- discretize(problem$elasticity)
  V <- disc$space; n <- V$ndof
  M1 <- problem$density * assemble_matrix(V, "mass")
  Z <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(n, n))
  M <- rbind(cbind(M1, Z), cbind(Z, M1))
  vec0 <- function(v) {
    if (is.numeric(v) && length(v) == 2L) return(c(rep(v[1], n), rep(v[2], n)))
    g <- eval_vector(v, V$dof_coords[, 1], V$dof_coords[, 2], "initial displacement")
    c(g[, 1], g[, 2])
  }
  vals <- disc$dirichlet$values
  out <- newmark(M, problem$damping * M, disc$A, function(t) problem$load_factor(t) * disc$b,
                 disc$dirichlet$dofs, function(t) vals, vec0(problem$initial), vec0(problem$velocity),
                 problem$time, dt, output_times, solver, ...)
  out$diagnostics$assembly_time <- disc$assembly_time
  split <- function(U) list(ux = U[seq_len(n), , drop = FALSE], uy = U[n + seq_len(n), , drop = FALSE])
  structure(list(times = out$times, values = split(out$U), velocity = split(out$Vel), space = V,
                 problem = problem, diagnostics = out$diagnostics), class = "pde_transient")
}

# Newmark average acceleration for M a + C v + K u = F(t) with prescribed
# values g(t) on `bdofs`. Boundary accelerations and velocities follow from
# the Newmark relations applied to the prescribed displacements.
newmark <- function(M, C, K, load, bdofs, g, u, v, time, dt, output_times, solver, ...) {
  check_number(dt, "dt", positive = TRUE)
  beta <- 0.25; gam <- 0.5
  t0 <- time[1]; t1 <- time[2]
  N <- max(1L, as.integer(ceiling((t1 - t0) / dt - 1e-10)))
  dt <- (t1 - t0) / N
  grid <- t0 + (0:N) * dt
  out_idx <- output_index(output_times, grid, dt)
  n <- nrow(M)
  free <- setdiff(seq_len(n), bdofs)
  if (length(bdofs)) u[bdofs] <- g(t0)
  ts <- proc.time()[["elapsed"]]
  max_rel <- 0; bad <- NULL
  note <- function(s) { max_rel <<- max(max_rel, s$info$relative_residual); if (!solve_ok(s$info)) bad <<- s$info }
  # Consistent initial acceleration on the free rows; on prescribed rows the
  # one-sided second difference of the data.
  a <- numeric(n)
  if (length(bdofs)) a[bdofs] <- (g(t0 + 2 * dt) - 2 * g(t0 + dt) + g(t0)) / dt^2
  r0 <- load(t0) - as.numeric(C %*% v) - as.numeric(K %*% u) - as.numeric(M[, bdofs, drop = FALSE] %*% a[bdofs])
  sM <- linear_solver(M[free, free, drop = FALSE], solver, ...)
  s0 <- sM$solve(r0[free]); note(s0)
  a[free] <- s0$x
  c0 <- 1 / (beta * dt^2); c1 <- gam / (beta * dt)
  Keff <- c0 * M + c1 * C + K
  slv <- linear_solver(Keff[free, free, drop = FALSE], solver, ...)
  Kfd <- Keff[free, bdofs, drop = FALSE]
  energy <- function(u, v) 0.5 * sum(v * as.numeric(M %*% v)) + 0.5 * sum(u * as.numeric(K %*% u))
  U <- Vel <- matrix(NA_real_, n, length(out_idx))
  E <- numeric(N + 1L); E[1] <- energy(u, v)
  k_out <- 1L
  if (out_idx[1] == 1L) { U[, 1] <- u; Vel[, 1] <- v; k_out <- 2L }
  for (k in seq_len(N)) {
    tn <- grid[k + 1L]
    rhs <- load(tn) + as.numeric(M %*% (c0 * u + dt * c0 * v + (1 / (2 * beta) - 1) * a)) +
      as.numeric(C %*% (c1 * u + (gam / beta - 1) * v + dt * (gam / (2 * beta) - 1) * a))
    un <- numeric(n)
    r <- rhs[free]
    if (length(bdofs)) { un[bdofs] <- g(tn); r <- r - as.numeric(Kfd %*% un[bdofs]) }
    s <- slv$solve(r, x0 = u[free]); note(s)
    un[free] <- s$x
    an <- c0 * (un - u - dt * v) - (1 / (2 * beta) - 1) * a
    v <- v + dt * ((1 - gam) * a + gam * an)
    u <- un; a <- an
    E[k + 1L] <- energy(u, v)
    if (k_out <= length(out_idx) && out_idx[k_out] == k + 1L) { U[, k_out] <- u; Vel[, k_out] <- v; k_out <- k_out + 1L }
  }
  if (!is.null(bad)) warn_solver(bad, "linear solve of at least one time step")
  diag <- structure(list(ndof = n, free = length(free), steps = N, dt = dt, method = "newmark",
                         solver = slv$info, energy = E, max_relative_residual = max_rel,
                         solve_time = proc.time()[["elapsed"]] - ts),
                    class = "pde_diagnostics")
  list(times = grid[out_idx], U = U, Vel = Vel, diagnostics = diag)
}
