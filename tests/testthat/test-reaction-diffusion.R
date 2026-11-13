test_that("implicit Euler Newton converges quadratically", {
  p <- reaction_diffusion_problem(pde_mesh_rectangle(c(0, 10), c(0, 1), nx = 20, ny = 2), diffusion = 0.1,
                                  reaction = function(u, x, y, t) u * (1 - u),
                                  initial = list(u = function(x, y) as.numeric(x < 1)), time = c(0, 0.2))
  s <- solve_pde(p, dt = 0.05, method = "implicit-euler")
  expect_lte(max(s$diagnostics$newton_iterations), 4)
})

test_that("linear decay matches the analytical solution", {
  D <- 0.1; c0 <- 2
  ex <- function(x, y, t) exp(-(2 * pi^2 * D + c0) * t) * sin(pi * x) * sin(pi * y)
  prob <- reaction_diffusion_problem(pde_mesh_rectangle(nx = 16), species = "u", diffusion = D,
                                     reaction = function(u, x, y, t) -c0 * u,
                                     initial = list(u = u_sin), bcs = list(u = dirichlet_bc(0)),
                                     time = c(0, 0.2))
  for (m in c("imex-euler", "implicit-euler")) {
    s <- solve_pde(prob, dt = 0.005, method = m)
    e <- errornorm(solution_at(s, 0.2), function(x, y) ex(x, y, 0.2), "L2")
    expect_lt(as.numeric(e), 5e-3)
  }
})

test_that("logistic growth with zero flux follows the ODE solution", {
  r <- 3; u0 <- 0.1; Tend <- 1
  exact <- 1 / (1 + (1 / u0 - 1) * exp(-r * Tend))
  prob <- reaction_diffusion_problem(pde_mesh_rectangle(nx = 4), diffusion = 1,
                                     reaction = function(u, x, y, t) r * u * (1 - u),
                                     initial = list(u = u0), time = c(0, Tend))
  e <- vapply(c(0.02, 0.01), function(dt) {
    s <- solve_pde(prob, dt = dt, method = "implicit-euler")
    max(abs(s$values$u[, length(s$times)] - exact))
  }, numeric(1))
  expect_lt(e[2], 0.01)
  expect_equal(log(e[1] / e[2], 2), 1, tolerance = 0.15)
  s <- solve_pde(prob, dt = 0.01, method = "implicit-euler")
  expect_true(max(s$diagnostics$newton_iterations) <= 6)
})

test_that("coupled linear systems match the matrix exponential", {
  A <- rbind(c(-1, 0.5), c(0.3, -2))
  prob <- reaction_diffusion_problem(pde_mesh_rectangle(nx = 3), species = c("a", "b"),
                                     diffusion = list(1, 0.2),
                                     reaction = function(u, x, y, t) u %*% t(A),
                                     initial = list(a = 1, b = 2), time = c(0, 0.5))
  # exact: u(t) = expm(A t) u0
  ev <- eigen(A); expm <- ev$vectors %*% diag(exp(ev$values * 0.5)) %*% solve(ev$vectors)
  exact <- Re(expm %*% c(1, 2))
  s <- solve_pde(prob, dt = 0.001, method = "implicit-euler", output_times = c(0, 0.5))
  last <- length(s$times)
  expect_equal(mean(s$values$a[, last]), exact[1], tolerance = 5e-3)
  expect_equal(mean(s$values$b[, last]), exact[2], tolerance = 5e-3)
  s2 <- solve_pde(prob, dt = 0.001, method = "implicit-euler", output_times = c(0, 0.5))
  expect_identical(s$values$a, s2$values$a)
  expect_s3_class(solution_at(s, 0.5, "b"), "fe_function")
})

test_that("user Jacobians are used and agree with finite differences", {
  reac <- function(u, x, y, t) cbind(u[, 1] * u[, 2] - u[, 1], -u[, 1] * u[, 2])
  jac <- function(u, x, y, t) {
    J <- array(0, c(nrow(u), 2, 2))
    J[, 1, 1] <- u[, 2] - 1; J[, 1, 2] <- u[, 1]
    J[, 2, 1] <- -u[, 2]; J[, 2, 2] <- -u[, 1]
    J
  }
  base <- list(mesh = pde_mesh_rectangle(nx = 4), species = c("p", "q"), diffusion = 0.1,
               reaction = reac, initial = list(p = function(x, y) 0.5 + x, q = 1), time = c(0, 0.3))
  s1 <- solve_pde(do.call(reaction_diffusion_problem, base), dt = 0.05, method = "implicit-euler")
  s2 <- solve_pde(do.call(reaction_diffusion_problem, c(base, jacobian = jac)), dt = 0.05,
                  method = "implicit-euler")
  expect_equal(s1$values$p, s2$values$p, tolerance = 1e-7)
  calls <- 0
  jac_count <- function(u, x, y, t) { calls <<- calls + 1; jac(u, x, y, t) }
  solve_pde(do.call(reaction_diffusion_problem, c(base, jacobian = jac_count)), dt = 0.05,
            method = "implicit-euler")
  expect_gt(calls, 0)
})
