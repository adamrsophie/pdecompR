test_that("Lagrange bases satisfy the nodal property and partition of unity", {
  for (p in 1:2) {
    el <- fe_element(p)
    B <- el$basis(el$nodes[, 1], el$nodes[, 2])
    expect_equal(B, diag(el$ndof), tolerance = 1e-14)
    xi <- runif(20) / 2; eta <- runif(20) / 2
    expect_equal(rowSums(el$basis(xi, eta)), rep(1, 20), tolerance = 1e-14)
    G <- el$grad(xi, eta)
    expect_equal(rowSums(G$xi), rep(0, 20), tolerance = 1e-13)
    expect_equal(rowSums(G$eta), rep(0, 20), tolerance = 1e-13)
  }
  expect_error(fe_element(3), "degree 1 and 2")
})

test_that("basis gradients match finite differences", {
  for (p in 1:2) {
    el <- fe_element(p)
    xi <- c(0.2, 0.1, 0.3); eta <- c(0.3, 0.6, 0.1); h <- 1e-6
    G <- el$grad(xi, eta)
    fdx <- (el$basis(xi + h, eta) - el$basis(xi - h, eta)) / (2 * h)
    fdy <- (el$basis(xi, eta + h) - el$basis(xi, eta - h)) / (2 * h)
    expect_equal(G$xi, fdx, tolerance = 1e-8)
    expect_equal(G$eta, fdy, tolerance = 1e-8)
  }
})

test_that("P2 spaces number vertex and edge dofs", {
  m <- pde_mesh_rectangle(nx = 2)
  V <- fe_space(m, 2)
  expect_equal(V$ndof, nrow(m$nodes) + nrow(m$edges))
  expect_equal(ncol(V$cell_dofs), 6)
  u <- interpolate(function(x, y) x^2 + x * y - y^2, V)
  pts <- cbind(runif(30), runif(30))
  expect_equal(evaluate(u, pts), pts[, 1]^2 + pts[, 1] * pts[, 2] - pts[, 2]^2, tolerance = 1e-12)
  expect_length(boundary_dofs(V, "bottom"), 5)
  expect_error(fe_space(m, family = "RT"), "Lagrange")
})
