test_that("Dirichlet Laplacian eigenvalues converge at the expected orders", {
  ex <- pi^2 * c(2, 5, 5, 8)
  for (deg in 1:2) {
    e <- vapply(c(8, 16), function(n) {
      ev <- solve_eigen(elliptic_problem(pde_mesh_rectangle(nx = n), bcs = dirichlet_bc(0), degree = deg), n = 4)
      expect_true(all(ev$residuals <= 1e-10))
      max(abs(ev$values - ex) / ex)
    }, numeric(1))
    expect_gt(log2(e[1] / e[2]), 2 * deg - 0.2)
  }
})

test_that("Neumann and elasticity eigenproblems", {
  ev <- solve_eigen(elliptic_problem(pde_mesh_rectangle(nx = 12)), n = 4)
  expect_lt(abs(ev$values[1]), 1e-10)
  expect_equal(ev$values[2:3] / pi^2, c(1, 1), tolerance = 0.01)
  V <- ev$vectors[[2]]$space; x <- ev$vectors[[2]]$values
  expect_equal(sum(x * as.numeric(assemble_matrix(V, "mass") %*% x)), 1, tolerance = 1e-10)
  ee <- solve_eigen(elasticity_problem(pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 8, ny = 2), E = 100,
                                       bcs = displacement_bc(0, "left")), n = 2)
  expect_true(all(ee$values > 0))
  expect_named(ee$vectors[[1]], c("ux", "uy"))
  set.seed(3); a <- stats::runif(1); set.seed(3)
  solve_eigen(elliptic_problem(pde_mesh_rectangle(nx = 4), bcs = dirichlet_bc(0)), n = 1)
  expect_equal(stats::runif(1), a)
  expect_error(solve_eigen(elliptic_problem(pde_mesh_rectangle(nx = 4), advection = c(1, 0))), "symmetric")
})
