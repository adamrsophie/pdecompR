# Singular solution on the L-shaped domain (-1, 1)^2 minus the quadrant
# x < 0, y > 0: u = r^(2/3) sin(2 phi / 3), phi measured from the edge on the
# negative x axis. It is harmonic and vanishes on the two reentrant edges.
l_vertices <- rbind(c(-1, -1), c(1, -1), c(1, 1), c(0, 1), c(0, 0), c(-1, 0))
l_phi <- function(x, y) (atan2(y, x) - pi) %% (2 * pi)
u_l <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  r^(2 / 3) * sin(2 * l_phi(x, y) / 3)
}
g_l <- function(x, y) {
  r <- sqrt(x^2 + y^2); th <- atan2(y, x); ph <- l_phi(x, y)
  a <- 2 / 3 * r^(-1 / 3)
  cbind(a * (sin(2 * ph / 3) * cos(th) - cos(2 * ph / 3) * sin(th)),
        a * (sin(2 * ph / 3) * sin(th) + cos(2 * ph / 3) * cos(th)))
}

test_that("the residual estimator is equivalent to the error under uniform refinement", {
  prob <- poisson_problem(pde_mesh_rectangle(nx = 4), source = f_sin, boundary = dirichlet_bc(0))
  est <- vapply(c(8, 16, 32), function(n) {
    s <- solve_pde(update_mesh(prob, pde_mesh_rectangle(nx = n)))
    c(estimate_error(s)$total, errornorm(s, u_sin, "H1_semi", gradient = g_sin))
  }, numeric(2))
  ratio <- est[1, ] / est[2, ]
  expect_true(all(ratio > 1 & ratio < 20))
  expect_lt(max(ratio) / min(ratio), 1.3)
  expect_equal(log2(est[1, 2] / est[1, 3]), 1, tolerance = 0.1)
  s <- solve_pde(update_mesh(prob, pde_mesh_rectangle(nx = 8)))
  e <- estimate_error(s)
  expect_length(e$eta, 128)
  expect_equal(e$total, sqrt(sum(e$eta^2)))
})

test_that("the estimator also works for P2, Neumann and variable diffusion", {
  prob <- elliptic_problem(pde_mesh_rectangle(nx = 8), diffusion = function(x, y) 1 + x, source = 1,
                           bcs = list(dirichlet_bc(0, "left"), neumann_bc(1, "right")), degree = 2)
  e <- estimate_error(solve_pde(prob))
  expect_true(all(is.finite(e$eta)) && e$total > 0)
})

test_that("Doerfler marking selects a minimal bulk set", {
  eta <- c(0.1, 0.5, 0.2, 0.05, 0.4)
  expect_equal(mark_elements(eta, 0.5), 2L)
  mk <- mark_elements(eta, 0.8)
  expect_true(sum(eta[mk]^2) >= 0.8 * sum(eta^2))
  expect_equal(mark_elements(eta, 0.7, "maximum"), c(2L, 5L))
  expect_error(mark_elements(eta, 1.5), "theta")
})

test_that("adaptive refinement beats uniform refinement on the L-shaped domain", {
  m0 <- pde_mesh_polygon(l_vertices, h = 0.5)
  # u vanishes on the two reentrant edges
  expect_lt(max(abs(u_l(c(-0.5, -0.2), c(0, 0)))), 1e-14)
  expect_lt(max(abs(u_l(c(0, 0), c(0.3, 0.9)))), 1e-12)
  prob <- poisson_problem(m0, source = 0, boundary = dirichlet_bc(u_l))
  ad <- adaptive_solve(prob, theta = 0.5, max_iterations = 12, max_dofs = 2500,
                       exact = u_l, gradient = g_l)
  h <- ad$history
  expect_true(all(diff(h$estimator) < 0))
  # adaptive rate in terms of ndof is close to the optimal -1/2
  k <- nrow(h)
  rate_ad <- log(h$error_H1_semi[k] / h$error_H1_semi[k - 3]) / log(h$ndof[k] / h$ndof[k - 3])
  expect_lt(rate_ad, -0.4)
  # uniform refinement with at least as many dofs has a larger error
  mu <- m0
  repeat {
    su <- solve_pde(update_mesh(prob, mu))
    if (su$space$ndof >= h$ndof[k]) break
    mu <- refine_uniform(mu)
  }
  eu <- as.numeric(errornorm(su, u_l, "H1_semi", gradient = g_l))
  expect_lt(h$error_H1_semi[k], eu)
  expect_s3_class(ad, "pde_adaptive")
  expect_output(print(ad), "Adaptive")
})
