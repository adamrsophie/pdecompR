test_that("advection-diffusion converges for moderate Peclet numbers", {
  b <- c(1, 0.5); k <- 0.1
  f <- function(x, y) {
    u <- u_sin(x, y); g <- g_sin(x, y)
    k * 2 * pi^2 * u + b[1] * g[, 1] + b[2] * g[, 2]
  }
  prob <- elliptic_problem(pde_mesh_rectangle(nx = 4), diffusion = k, advection = b, source = f,
                           bcs = dirichlet_bc(0))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.8)
  sol <- solve_pde(prob)
  expect_equal(sol$diagnostics$solver$method, "lu")
  # SUPG is consistent and also converges
  ps <- elliptic_problem(prob$mesh, diffusion = k, advection = b, source = f, bcs = dirichlet_bc(0),
                         stabilization = "supg")
  cs2 <- convergence_study(ps, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs2$rate_H1_semi[-1]), 0.9)
})

test_that("SUPG removes the oscillations of Galerkin at high Peclet number", {
  eps <- 1e-3
  m <- pde_mesh_rectangle(nx = 16)
  common <- list(mesh = m, diffusion = eps, advection = c(1, 0), source = 1, bcs = dirichlet_bc(0))
  g <- solve_pde(do.call(elliptic_problem, common))
  s <- solve_pde(do.call(elliptic_problem, c(common, stabilization = "supg")))
  # exact solution is bounded by 1 (u ~ x away from the outflow layer)
  over_g <- max(g$u$values) - 1
  over_s <- max(s$u$values) - 1
  expect_gt(over_g, 0.1)
  expect_lt(over_s, 0.05)
  expect_lt(min(g$u$values), -0.1)
  expect_gt(min(s$u$values), -0.05)
  expect_error(elliptic_problem(m, advection = c(1, 0), stabilization = "supg", degree = 2), "P1")
})

test_that("transient advection-diffusion with SUPG is consistent", {
  b <- c(1, 0); k <- 0.05
  ex <- function(x, y, t) exp(-t) * sin(pi * x) * sin(pi * y)
  src <- function(x, y, t) {
    u <- ex(x, y, t)
    -u + k * 2 * pi^2 * u + b[1] * pi * exp(-t) * cos(pi * x) * sin(pi * y)
  }
  e <- vapply(c("none", "supg"), function(st) {
    prob <- parabolic_problem(pde_mesh_rectangle(nx = 16), diffusion = k, advection = b, source = src,
                              initial = u_sin, bcs = dirichlet_bc(0), time = c(0, 0.5),
                              stabilization = st)
    s <- solve_pde(prob, dt = 0.01, output_times = 0.5)
    as.numeric(errornorm(solution_at(s, 0.5), function(x, y) ex(x, y, 0.5), "L2"))
  }, numeric(1))
  expect_true(all(e < 5e-3))
})
