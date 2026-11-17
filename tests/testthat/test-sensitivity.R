test_that("forward sensitivities match the analytic derivative", {
  m <- pde_mesh_rectangle(nx = 8)
  pf <- function(theta) poisson_problem(m, diffusion = theta[1], source = 1, boundary = dirichlet_bc(0))
  s <- pde_sensitivity(pf, 2)
  expect_equal(s$sensitivities[, 1], -s$solution$u$values / 2, tolerance = 1e-8)
})

test_that("forward, adjoint and finite-difference gradients agree", {
  m <- pde_mesh_rectangle(nx = 8)
  pf <- function(theta) elliptic_problem(m, diffusion = function(x, y) theta[1] + theta[2] * x,
                                         reaction = theta[3],
                                         source = function(x, y) 1 + theta[2] * y,
                                         bcs = list(dirichlet_bc(function(x, y) theta[3] * (1 + y), "left"),
                                                    dirichlet_bc(0, "right")))
  V <- fe_space(m)
  obs <- pde_observations(c(0, 0.3, 0.6, 0.5), c(0.5, 0.4, 0.5, 0.8), c(0.4, 0.05, 0.1, 0.02), sd = 0.1)
  J <- misfit_objective(V, obs)
  th <- c(1, 0.5, 0.3)
  sf <- pde_sensitivity(pf, th, objective = J, method = "forward")
  sa <- pde_sensitivity(pf, th, objective = J, method = "adjoint")
  expect_equal(sf$gradient, sa$gradient, tolerance = 1e-8)
  fd <- vapply(1:3, function(j) {
    h <- 1e-5; tp <- th; tm <- th; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (J(solve_pde(pf(tp))$u$values)$value - J(solve_pde(pf(tm))$u$values)$value) / (2 * h)
  }, numeric(1))
  expect_equal(sa$gradient, fd, tolerance = 1e-6)
  expect_error(pde_sensitivity(pf, th, method = "adjoint"), "objective")
})

test_that("misfit objectives have consistent gradients", {
  V <- fe_space(pde_mesh_rectangle(nx = 4))
  obs <- pde_observations(c(0.2, 0.7), c(0.3, 0.9), c(1, 2), sd = c(0.5, 1))
  J <- misfit_objective(V, obs)
  u <- runif(V$ndof)
  g <- J(u)$gradient
  e <- numeric(V$ndof); e[7] <- 1e-6
  expect_equal((J(u + e)$value - J(u - e)$value) / 2e-6, g[7], tolerance = 1e-6)
  expect_error(pde_observations(1, 1:2, 1), "equal length")
  expect_error(pde_observations(1, 1, 1, sd = 0), "positive")
})

test_that("difference steps respect parameter bounds", {
  m <- pde_mesh_rectangle(nx = 4)
  pf <- function(theta) {
    if (theta[1] <= 0) stop("diffusion must be positive")
    poisson_problem(m, diffusion = theta[1], source = 1, boundary = dirichlet_bc(0))
  }
  s <- pde_sensitivity(pf, 1e-8, lower = 1e-8)
  expect_equal(s$sensitivities[, 1], -s$solution$u$values / 1e-8, tolerance = 1e-4)
})

test_that("nonlinear and transient adjoint gradients pass Taylor and difference checks", {
  mesh <- pde_mesh_rectangle(nx = 5)
  ob <- misfit_objective(fe_space(mesh), pde_observations(c(0.5, 0.3), c(0.5, 0.6), c(0.05, 0.02)))
  pf <- function(th) nonlinear_problem(mesh, diffusion = function(u) th[1] + u^2, source = function(x, y) th[2] * (1 + x),
                                       reaction = function(u) th[3] * u^3, bcs = dirichlet_bc(function(x, y) th[3] * x * 0.1))
  J <- function(th) ob(solve_pde(pf(th), tol = 1e-13)$u$values)$value
  th <- c(1, 2, 0.5)
  s <- nonlinear_sensitivity(pf, th, ob, tol = 1e-13)
  d <- c(1, -1, 0.5)
  rem <- vapply(10^-(2:3), function(e) abs(J(th + e * d) - s$value - e * sum(s$gradient * d)), numeric(1))
  expect_gt(log10(rem[1] / rem[2]), 1.9)
  pt <- function(th) parabolic_problem(mesh, diffusion = th[1], reaction = th[2],
                                       source = function(x, y, t) th[3] * sin(t) * x,
                                       initial = function(x, y) th[4] * sin(pi * x) * sin(pi * y),
                                       bcs = dirichlet_bc(function(x, y, t) 0.1 * t * y), time = c(0, 0.3))
  for (me in c("backward-euler", "crank-nicolson")) {
    Jt <- function(th) { s <- solve_pde(pt(th), dt = 0.05, method = me, startup = 0L); ob(s$values[, ncol(s$values)])$value }
    th <- c(1, 0.5, 2, 1)
    g <- transient_sensitivity(pt, th, ob, dt = 0.05, method = me)$gradient
    fd <- vapply(1:4, function(j) { e <- replace(numeric(4), j, 1e-5); (Jt(th + e) - Jt(th - e)) / 2e-5 }, numeric(1))
    expect_equal(g, fd, tolerance = 1e-6)
  }
  expect_error(transient_sensitivity(function(th) poisson_problem(mesh, source = th), 1, ob, dt = 0.1), "parabolic")
})
