exact_monomial <- function(i, j) factorial(i) * factorial(j) / factorial(i + j + 2)

test_that("triangle rules integrate polynomials of their degree exactly", {
  for (d in 0:10) {
    q <- quadrature_triangle(d)
    expect_gte(q$degree, d)
    expect_equal(sum(q$weights), 0.5, tolerance = 1e-14)
    expect_true(all(q$points >= 0) && all(rowSums(q$points) <= 1))
    for (i in 0:d) for (j in 0:(d - i)) {
      expect_equal(sum(q$weights * q$points[, 1]^i * q$points[, 2]^j), exact_monomial(i, j),
                   tolerance = 1e-13, label = sprintf("degree %d monomial x^%d y^%d", d, i, j))
    }
  }
})

test_that("rules are not exact beyond their degree", {
  q <- quadrature_triangle(2)
  expect_false(isTRUE(all.equal(sum(q$weights * q$points[, 1]^3), exact_monomial(3, 0))))
})

test_that("Gauss-Legendre rules are exact to degree 2n - 1", {
  for (n in 1:8) {
    q <- quadrature_interval(n)
    for (k in 0:(2 * n - 1)) {
      expect_equal(sum(q$weights * q$points[, 1]^k), 1 / (k + 1), tolerance = 1e-13)
    }
  }
})

test_that("integration over a mesh maps to physical elements", {
  m <- pde_mesh_rectangle(c(1, 3), c(-1, 2), nx = 3, ny = 2, pattern = "alternating")
  expect_equal(integrate_mesh(m, function(x, y) x^2 * y), (27 - 1) / 3 * (4 - 1) / 2, tolerance = 1e-12)
  L <- pde_mesh_polygon(rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2)))
  expect_equal(integrate_mesh(L, function(x, y) 1), 3)
})
