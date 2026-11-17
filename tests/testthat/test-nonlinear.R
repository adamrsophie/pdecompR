test_that("quasilinear diffusion converges with Newton's method", {
  k <- function(u) 1 + u^2
  f <- function(x, y) {
    u <- u_sin(x, y); g <- g_sin(x, y)
    k(u) * 2 * pi^2 * u - 2 * u * (g[, 1]^2 + g[, 2]^2)
  }
  prob <- nonlinear_problem(pde_mesh_rectangle(nx = 4), diffusion = k, diffusion_derivative = function(u) 2 * u,
                            source = f, bcs = dirichlet_bc(0))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.9)
  expect_gt(min(cs$rate_H1_semi[-1]), 0.95)
  s <- solve_pde(update_mesh(prob, pde_mesh_rectangle(nx = 16)))
  expect_true(s$converged)
  r <- s$newton$residual
  # quadratic convergence: each step roughly squares the relative residual
  rr <- r / r[1]
  idx <- which(rr < 1e-2 & rr > 1e-13)
  if (length(idx) >= 2L) expect_gt(min(log(rr[idx[-1]]) / log(rr[idx[-length(idx)]])), 1.6)
  expect_lt(tail(r, 1), 1e-10 * r[1])
  expect_lte(nrow(s$newton), 8)
  # finite-difference derivative gives the same solution
  prob2 <- prob; prob2$diffusion_derivative <- NULL
  s2 <- solve_pde(update_mesh(prob2, pde_mesh_rectangle(nx = 16)))
  expect_equal(s$u$values, s2$u$values, tolerance = 1e-8)
})

test_that("semilinear reaction (Bratu-type) problems converge", {
  # -Laplace u + u^3 = f with u = sin sin
  f <- function(x, y) f_sin(x, y) + u_sin(x, y)^3
  prob <- nonlinear_problem(pde_mesh_rectangle(nx = 4), reaction = function(u) u^3,
                            reaction_derivative = function(u) 3 * u^2, source = f, bcs = dirichlet_bc(0))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.9)
})

test_that("regularized p-Laplacian converges to a manufactured solution", {
  p <- 3; eps <- 0.1
  # f = -div((|grad u|^2 + eps^2)^((p-2)/2) grad u) for u = sin sin
  f <- function(x, y) {
    ux <- pi * cos(pi * x) * sin(pi * y); uy <- pi * sin(pi * x) * cos(pi * y)
    uxx <- -pi^2 * sin(pi * x) * sin(pi * y); uyy <- uxx; uxy <- pi^2 * cos(pi * x) * cos(pi * y)
    s <- ux^2 + uy^2 + eps^2
    k <- s^((p - 2) / 2)
    dk <- (p - 2) * s^((p - 4) / 2)
    -(k * (uxx + uyy) + dk * (ux^2 * uxx + 2 * ux * uy * uxy + uy^2 * uyy))
  }
  prob <- p_laplace_problem(pde_mesh_rectangle(nx = 4), p = p, epsilon = eps, source = f,
                            bcs = dirichlet_bc(0))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.85)
  expect_gt(min(cs$rate_H1_semi[-1]), 0.95)
  pl2 <- p_laplace_problem(pde_mesh_rectangle(nx = 8), p = 2, source = f_sin, bcs = dirichlet_bc(0))
  s2 <- solve_pde(pl2)
  s0 <- solve_pde(poisson_problem(pl2$mesh, source = f_sin, boundary = dirichlet_bc(0)))
  expect_equal(s2$u$values, s0$u$values, tolerance = 1e-10)
  expect_error(p_laplace_problem(pl2$mesh, p = 1), "larger than 1")
})

test_that("Newton failure is reported, not hidden", {
  prob <- nonlinear_problem(pde_mesh_rectangle(nx = 4), diffusion = function(u) 1 + u^2,
                            source = 50, bcs = dirichlet_bc(0))
  expect_warning(s <- solve_pde(prob, maxit = 1), "Newton stopped")
  expect_false(s$converged)
})
