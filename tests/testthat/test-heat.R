test_that("Crank-Nicolson and backward Euler have temporal orders 2 and 1", {
  # u = exp(-t) (x + 2 y) is P1 exact in space, so only time error remains.
  ex <- function(x, y, t) exp(-t) * (x + 2 * y)
  m <- pde_mesh_rectangle(nx = 4)
  prob <- heat_problem(m, source = function(x, y, t) -exp(-t) * (x + 2 * y),
                       initial = function(x, y) x + 2 * y, bcs = dirichlet_bc(ex), time = c(0, 1))
  err <- function(dt, method) {
    s <- solve_pde(prob, dt = dt, method = method, output_times = 1)
    xy <- s$space$dof_coords
    max(abs(s$values[, 1] - ex(xy[, 1], xy[, 2], 1)))
  }
  dts <- c(0.05, 0.025, 0.0125)
  e_cn <- vapply(dts, err, numeric(1), method = "crank-nicolson")
  e_be <- vapply(dts, err, numeric(1), method = "backward-euler")
  r_cn <- diff(log(e_cn)) / diff(log(dts))
  r_be <- diff(log(e_be)) / diff(log(dts))
  expect_true(all(r_cn > 1.9 & r_cn < 2.1))
  expect_true(all(r_be > 0.9 & r_be < 1.1))
})

test_that("heat equation decays like the analytical solution", {
  kappa <- 0.5
  ex <- function(x, y, t) exp(-2 * pi^2 * kappa * t) * sin(pi * x) * sin(pi * y)
  prob <- heat_problem(pde_mesh_rectangle(nx = 24), diffusivity = kappa, initial = u_sin,
                       bcs = dirichlet_bc(0), time = c(0, 0.1))
  s <- solve_pde(prob, dt = 0.0025)
  expect_equal(length(s$times), 41)
  u_end <- solution_at(s, 0.1)
  e <- errornorm(u_end, function(x, y) ex(x, y, 0.1), "L2")
  expect_lt(as.numeric(e), 2e-3)
  expect_output(print(s), "pde_transient")
  expect_error(solve_pde(prob, dt = 0.01, output_times = 0.005), "time grid")
})

test_that("time-dependent Neumann data and sources are handled", {
  # u = t * (x^2) has du/dt = x^2, -Laplace u = -2 t, flux on right = 2 t
  ex <- function(x, y, t) t * x^2
  prob <- heat_problem(pde_mesh_rectangle(nx = 8), source = function(x, y, t) x^2 - 2 * t,
                       initial = 0, bcs = list(dirichlet_bc(ex, c("left", "top", "bottom")),
                                               neumann_bc(function(x, y, t) 2 * t, "right")),
                       time = c(0, 0.5))
  s <- solve_pde(prob, dt = 0.05, output_times = c(0.25, 0.5))
  e <- errornorm(solution_at(s, 0.5), function(x, y) ex(x, y, 0.5), "L2")
  expect_lt(as.numeric(e), 5e-3)
})

test_that("mass matrix conservation holds for zero-flux heat problems", {
  prob <- heat_problem(pde_mesh_rectangle(nx = 8), initial = function(x, y) exp(-20 * ((x - 0.3)^2 + y^2)),
                       time = c(0, 0.2))
  s <- solve_pde(prob, dt = 0.02, method = "backward-euler")
  w <- assemble_vector(s$space, 1)
  tot <- colSums(w * s$values)
  expect_lt(max(abs(tot - tot[1])), 1e-12)
  # backward Euler is monotone in time for the energy
  M <- assemble_matrix(s$space, "mass")
  en <- apply(s$values, 2, function(u) sum(u * (M %*% u)))
  expect_true(all(diff(en) <= 1e-14))
})

test_that("BDF2 and SDIRK2 are second order and adaptive SDIRK2 meets its tolerance", {
  m <- pde_mesh_rectangle(nx = 6)
  p <- heat_problem(m, initial = function(x, y) sin(pi * x) * sin(pi * y),
                    source = function(x, y, t) cos(3 * t) * x * (1 - x), bcs = dirichlet_bc(0), time = c(0, 1))
  ref <- solution_at(solve_pde(p, dt = 1 / 1024, method = "sdirk2", output_times = 1), 1)
  for (me in c("bdf2", "sdirk2")) {
    e <- vapply(c(32, 64), function(N) errornorm(solution_at(solve_pde(p, dt = 1 / N, method = me, output_times = 1), 1), ref),
                numeric(1))
    expect_gt(log2(e[1] / e[2]), 1.85)
  }
  s <- solve_pde(p, dt = 0.1, method = "sdirk2", adaptive = TRUE, tol = 1e-5, output_times = c(0.25, 1))
  expect_equal(s$times, c(0.25, 1))
  expect_lt(errornorm(solution_at(s, 1), ref), 5e-6)
  expect_gt(s$diagnostics$steps, 20)
  expect_error(solve_pde(p, dt = 0.1, adaptive = TRUE), "sdirk2")
  pt <- parabolic_problem(m, diffusion = function(x, y, t) 1 + t, initial = 0, bcs = dirichlet_bc(0), time = c(0, 1))
  expect_error(solve_pde(pt, dt = 0.1, method = "bdf2"), "do not depend on time")
})
