test_that("parameter estimation recovers diffusion coefficients", {
  m <- pde_mesh_rectangle(nx = 8)
  pf <- function(theta, mesh) elliptic_problem(mesh, diffusion = function(x, y) theta[1] + theta[2] * x,
                                               source = 1, bcs = dirichlet_bc(0))
  truth <- c(0.5, 1.5)
  sol <- solve_pde(pf(truth, m))
  set.seed(42)
  xo <- runif(20, 0.1, 0.9); yo <- runif(20, 0.1, 0.9)
  obs <- pde_observations(xo, yo, evaluate(sol$u, xo, yo), sd = 1e-3)
  est <- estimate_parameters(pf, c(1, 1), obs, mesh = m, lower = c(0.05, 0))
  expect_equal(est$estimate, truth, tolerance = 1e-4)
  expect_equal(est$convergence, 0L)
  expect_true(all(is.finite(est$std_error)))
  expect_false(is.null(est$discretization_check))
  expect_output(print(est), "measurement noise only")
  # regularization pulls the estimate toward the reference
  est_r <- estimate_parameters(pf, c(1, 1), obs, mesh = m, lower = c(0.05, 0), regularization = 1e6,
                               theta_ref = c(1, 1))
  expect_lt(sum((est_r$estimate - c(1, 1))^2), sum((est$estimate - c(1, 1))^2))
})

test_that("unconstrained optimal control matches the reduced gradient optimum", {
  m <- pde_mesh_rectangle(nx = 8)
  tgt <- function(x, y) sin(pi * x) * sin(pi * y)
  a <- optimal_control(m, tgt, alpha = 1e-3)
  expect_lt(a$gradient_norm, 1e-10)
  b <- optimal_control(m, tgt, alpha = 1e-3, lower = -1e6, upper = 1e6, control = list(maxit = 2000, factr = 1e2, pgtol = 0))
  expect_equal(b$objective, a$objective, tolerance = 1e-5)
  expect_equal(a$control$values, b$control$values, tolerance = 1e-2)
  # for tiny alpha the state approaches the target up to the P1 interpolation error
  c1 <- optimal_control(m, tgt, alpha = 1e-8)
  expect_lt(as.numeric(errornorm(c1$state, tgt, "L2")), 0.03)
  expect_lt(as.numeric(errornorm(c1$state, tgt, "L2")), as.numeric(errornorm(a$state, tgt, "L2")))
})

test_that("bounded control respects the bounds", {
  m <- pde_mesh_rectangle(nx = 8)
  oc <- optimal_control(m, function(x, y) 2 * sin(pi * x) * sin(pi * y), alpha = 1e-4, lower = 0, upper = 10)
  expect_true(all(oc$control$values >= -1e-12 & oc$control$values <= 10 + 1e-12))
  expect_true(any(oc$control$values >= 10 - 1e-8))
  expect_output(print(oc), "optimal control")
})
