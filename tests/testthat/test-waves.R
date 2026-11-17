test_that("Newmark wave solutions converge and conserve energy", {
  u0 <- function(x, y) sin(pi * x) * sin(pi * y); w <- pi * sqrt(2)
  e <- vapply(c(8, 16), function(n) {
    s <- solve_pde(wave_problem(pde_mesh_rectangle(nx = n), initial = u0, bcs = dirichlet_bc(0), time = c(0, 1)),
                   dt = 1 / (4 * n), output_times = 1)
    E <- s$diagnostics$energy
    expect_lt(diff(range(E)) / E[1], 1e-11)
    errornorm(solution_at(s, 1), function(x, y) cos(w) * u0(x, y))
  }, numeric(1))
  expect_gt(log2(e[1] / e[2]), 1.9)
  m <- pde_mesh_rectangle(nx = 6)
  p <- wave_problem(m, initial = u0, velocity = function(x, y) x * (1 - x) * y, source = function(x, y, t) sin(2 * t),
                    bcs = dirichlet_bc(function(x, y, t) 0.1 * sin(t) * x), time = c(0, 1))
  ref <- solution_at(solve_pde(p, dt = 1 / 1024, output_times = 1), 1)
  et <- vapply(c(32, 64), function(N) errornorm(solution_at(solve_pde(p, dt = 1 / N, output_times = 1), 1), ref), numeric(1))
  expect_gt(log2(et[1] / et[2]), 1.9)
  E <- solve_pde(wave_problem(m, initial = u0, damping = 0.5, bcs = dirichlet_bc(0), time = c(0, 1)), dt = 0.02)$diagnostics$energy
  expect_true(all(diff(E) <= 1e-12))
  expect_error(wave_problem(m, damping = -1), "nonnegative")
})

test_that("elastodynamics conserves energy and damps to the static solution", {
  el <- elasticity_problem(pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 8, ny = 2), E = 100, nu = 0.3,
                           body_force = c(0, -1), bcs = displacement_bc(0, "left"))
  s <- solve_pde(elastodynamics_problem(el, initial = function(x, y) cbind(0, -0.01 * x^2), time = c(0, 1),
                                        load_factor = function(t) 0), dt = 0.01)
  E <- s$diagnostics$energy
  expect_lt(diff(range(E)) / E[1], 1e-11)
  st <- solve_pde(el)
  d <- solve_pde(elastodynamics_problem(el, damping = 1, time = c(0, 40)), dt = 0.05, output_times = 40)
  expect_lt(max(abs(solution_at(d, 40, "uy")$values - st$uy$values)) / max(abs(st$uy$values)), 1e-5)
})
