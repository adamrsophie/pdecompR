test_that("quadrature assembly matches closed-form P1 element matrices", {
  P <- rbind(c(0.1, 0.2), c(1.3, 0.4), c(0.5, 1.7))
  m <- pde_mesh(P, rbind(1:3))
  V <- fe_space(m)
  K <- as.matrix(assemble_matrix(V, "stiffness", 2.5))
  M <- as.matrix(assemble_matrix(V, "mass"))
  expect_equal(K, element_stiffness_p1(P, 2.5), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(M, element_mass_p1(P), tolerance = 1e-12, ignore_attr = TRUE)
  area <- mesh_quality(m)$area
  expect_equal(element_mass_p1(P), area / 12 * matrix(c(2, 1, 1, 1, 2, 1, 1, 1, 2), 3))
})

test_that("global matrices are sparse, symmetric and have the right kernel", {
  V <- fe_space(pde_mesh_rectangle(nx = 6))
  K <- assemble_matrix(V, "stiffness")
  M <- assemble_matrix(V, "mass")
  expect_s4_class(K, "dgCMatrix")
  expect_true(Matrix::isSymmetric(K, tol = 1e-13))
  expect_true(Matrix::isSymmetric(M, tol = 1e-13))
  expect_equal(max(abs(K %*% rep(1, V$ndof))), 0, tolerance = 1e-12)
  expect_equal(sum(M), 1, tolerance = 1e-13)
  ev <- eigen(as.matrix(M), only.values = TRUE)$values
  expect_gt(min(ev), 0)
  evk <- eigen(as.matrix(K), only.values = TRUE)$values
  expect_gt(min(evk), -1e-12)
  expect_equal(sum(abs(evk) < 1e-10), 1)
  # fill: at most 7 nonzeros per row for this structured P1 mesh
  expect_lte(max(diff(K@p)), 7)
  # energy of a linear function: int |grad u|^2 = |(2, -1)|^2 * area
  u <- interpolate(function(x, y) 2 * x - y, V)$values
  expect_equal(sum(u * as.numeric(K %*% u)), 5, tolerance = 1e-12)
  expect_equal(sum(u * as.numeric(M %*% u)), integrate_mesh(V$mesh, function(x, y) (2 * x - y)^2), tolerance = 1e-12)
})

test_that("assembly is invariant under node renumbering", {
  m <- pde_mesh_rectangle(nx = 4, pattern = "alternating")
  set.seed(1)
  perm <- sample(nrow(m$nodes))
  inv <- order(perm)
  m2 <- pde_mesh(m$nodes[perm, ], matrix(inv[m$triangles], ncol = 3),
                 cbind(inv[m$boundary$n1], inv[m$boundary$n2]), m$boundary$marker)
  K1 <- as.matrix(assemble_matrix(fe_space(m), "stiffness", function(x, y) 1 + x * y))
  K2 <- as.matrix(assemble_matrix(fe_space(m2), "stiffness", function(x, y) 1 + x * y))
  expect_equal(K2, K1[perm, perm], tolerance = 1e-13)
})

test_that("load vectors and advection matrices integrate correctly", {
  V <- fe_space(pde_mesh_rectangle(nx = 5))
  b <- assemble_vector(V, function(x, y) x + y)
  expect_equal(sum(b), 1, tolerance = 1e-13)
  C <- assemble_matrix(V, "advection", c(1, 0.5))
  u <- interpolate(function(x, y) x * y, V)$values
  # int (b . grad u) * 1 = int (y + 0.5 x) = 0.75
  expect_equal(sum(C %*% u), 0.75, tolerance = 1e-12)
  V2 <- fe_space(V$mesh, 2)
  expect_equal(sum(assemble_vector(V2, function(x, y) x^2)), 1 / 3, tolerance = 1e-13)
  expect_equal(sum(assemble_matrix(V2, "mass")), 1, tolerance = 1e-13)
})

test_that("boundary integrals use edge lengths and markers", {
  V <- fe_space(pde_mesh_rectangle(c(0, 2), c(0, 1), nx = 4, ny = 3))
  expect_equal(sum(assemble_boundary_vector(V, 1)), 6, tolerance = 1e-13)
  expect_equal(sum(assemble_boundary_vector(V, function(x, y) x, "bottom")), 2, tolerance = 1e-13)
  expect_equal(sum(assemble_boundary_matrix(V, 1, "left")), 1, tolerance = 1e-13)
  V2 <- fe_space(V$mesh, 2)
  expect_equal(sum(assemble_boundary_vector(V2, function(x, y) x^2, "top")), 8 / 3, tolerance = 1e-13)
})

test_that("tensor diffusion assembles grad v^T K grad u", {
  V <- fe_space(pde_mesh_rectangle(nx = 3))
  Kt <- matrix(c(2, 0.5, 0.5, 1), 2)
  A <- assemble_matrix(V, "stiffness", Kt)
  u <- interpolate(function(x, y) x + 2 * y, V)$values
  g <- c(1, 2)
  expect_equal(sum(u * as.numeric(A %*% u)), sum(g * (Kt %*% g)), tolerance = 1e-12)
  expect_error(assemble_matrix(V, "stiffness", matrix(c(1, 0, 0, -1), 2)), "positive definite")
})
