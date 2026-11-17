test_that("split parabolic runs reproduce one run for every stepper", {
  m <- pde_mesh_rectangle(nx = 4)
  mk <- function(T, td = FALSE) {
    parabolic_problem(m, diffusion = if (td) function(x, y, t) 1 + t * x else 1,
                      source = function(x, y, t) sin(3 * t) * x, initial = function(x, y) as.numeric(x > 0.5),
                      bcs = dirichlet_bc(function(x, y, t) t * y), time = c(0, T))
  }
  for (me in c("crank-nicolson", "backward-euler", "bdf2", "sdirk2")) {
    for (td in if (me %in% c("bdf2", "sdirk2")) FALSE else c(FALSE, TRUE)) {
      full <- solve_pde(mk(1, td), dt = 0.05, method = me, startup = 4)
      a <- solve_pde(mk(0.1, td), dt = 0.05, method = me, startup = 4)
      b <- solve_pde(mk(1, td), dt = 0.05, method = me, startup = 4, restart = a)
      expect_equal(b$state$u, full$state$u, tolerance = 1e-12)
      expect_equal(b$state$step, 20L)
    }
  }
  ad <- solve_pde(mk(1), dt = 0.05, method = "sdirk2", adaptive = TRUE, tol = 1e-4, output_times = c(0.5, 1))
  a1 <- solve_pde(mk(0.5), dt = 0.05, method = "sdirk2", adaptive = TRUE, tol = 1e-4)
  b1 <- solve_pde(mk(1), dt = 0.05, method = "sdirk2", adaptive = TRUE, tol = 1e-4, restart = a1)
  expect_equal(b1$state$u, ad$state$u, tolerance = 1e-12)
  f <- tempfile(fileext = ".rds"); saveRDS(a1, f)
  expect_equal(solve_pde(mk(1), dt = 0.05, method = "sdirk2", adaptive = TRUE, tol = 1e-4,
                         restart = readRDS(f))$state$u, b1$state$u)
  n <- 0; times <- numeric(0)
  invisible(solve_pde(mk(1), dt = 0.1, callback = function(t, u) { n <<- n + 1; times <<- c(times, t) }))
  expect_equal(n, 10); expect_equal(times[10], 1)
  ta <- numeric(0)
  ad2 <- solve_pde(mk(1), dt = 0.05, method = "sdirk2", adaptive = TRUE, tol = 1e-4,
                   callback = function(t, u) ta <<- c(ta, t))
  expect_length(ta, ad2$diagnostics$steps); expect_equal(ta[length(ta)], 1); expect_true(all(diff(ta) > 0))
  expect_error(solve_pde(mk(1), dt = 0.1, restart = solve_pde(mk(1), dt = 0.1)), "restart time")
  expect_error(solve_pde(mk(1), dt = 0.1, callback = 1), "callback")
})

test_that("split wave, elastodynamics and reaction-diffusion runs reproduce one run", {
  m <- pde_mesh_rectangle(nx = 4)
  wp <- function(T) wave_problem(m, damping = 0.3, source = function(x, y, t) cos(t) * x,
                                 initial = function(x, y) sin(pi * x) * sin(pi * y), velocity = function(x, y) x * (1 - x),
                                 bcs = dirichlet_bc(function(x, y, t) 0.1 * t * y), time = c(0, T))
  f <- solve_pde(wp(1), dt = 0.05); b <- solve_pde(wp(1), dt = 0.05, restart = solve_pde(wp(0.4), dt = 0.05))
  expect_equal(b$state$u, f$state$u, tolerance = 1e-12); expect_equal(b$state$v, f$state$v, tolerance = 1e-12)
  el <- elasticity_problem(pde_mesh_rectangle(c(0, 2), c(0, 1), nx = 4, ny = 2), E = 10, nu = 0.3, body_force = c(0, -1),
                           bcs = displacement_bc(c(0, 0), "left"))
  ep <- function(T) elastodynamics_problem(el, damping = 0.2, time = c(0, T), load_factor = function(t) sin(t))
  f <- solve_pde(ep(2), dt = 0.1); b <- solve_pde(ep(2), dt = 0.1, restart = solve_pde(ep(0.7), dt = 0.1))
  expect_equal(b$state$u, f$state$u, tolerance = 1e-12)
  rp <- function(T) reaction_diffusion_problem(m, c("a", "b"), diffusion = c(1, 0.5),
                                               reaction = function(U, x, y, t) cbind(-U[, 1] * U[, 2], U[, 1] * U[, 2]),
                                               initial = list(a = function(x, y) 1 + x, b = function(x, y) y), time = c(0, T))
  for (me in c("imex-euler", "implicit-euler")) {
    f <- solve_pde(rp(1), dt = 0.05, method = me)
    b <- solve_pde(rp(1), dt = 0.05, method = me, restart = solve_pde(rp(0.3), dt = 0.05, method = me))
    expect_equal(b$state$u, f$state$u, tolerance = 1e-12)
  }
  expect_error(solve_pde(wp(1), dt = 0.1, restart = b), "different kind")
})
