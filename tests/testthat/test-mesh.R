test_that("structured rectangle meshes have the expected size and markers", {
  for (pat in c("right", "left", "alternating")) {
    m <- pde_mesh_rectangle(c(0, 2), c(0, 1), nx = 4, ny = 3, pattern = pat)
    expect_equal(nrow(m$nodes), 5 * 4)
    expect_equal(nrow(m$triangles), 2 * 4 * 3)
    expect_equal(nrow(m$boundary), 2 * (4 + 3))
    expect_equal(sum(mesh_quality(m)$area), 2)
    expect_true(all(mesh_quality(m)$area > 0))
    expect_equal(sort(boundary_markers(m)), sort(c("bottom", "right", "top", "left")))
    # Euler characteristic of a disk: V - E + F = 1
    expect_equal(nrow(m$nodes) - nrow(m$edges) + nrow(m$triangles), 1)
  }
})

test_that("boundary nodes are selected by marker", {
  m <- pde_mesh_rectangle(nx = 3, ny = 2)
  left <- boundary_nodes(m, "left")
  expect_true(all(m$nodes[left, 1] == 0))
  expect_length(left, 3)
  expect_length(boundary_nodes(m), 2 * (3 + 2))
  expect_error(boundary_nodes(m, "inlet"), "Unknown boundary marker")
})

test_that("boundary edges are oriented with the domain on the left", {
  m <- pde_mesh_rectangle(nx = 3)
  P <- m$nodes; bd <- m$boundary
  mx <- (P[bd$n1, 1] + P[bd$n2, 1]) / 2; my <- (P[bd$n1, 2] + P[bd$n2, 2]) / 2
  nx <- P[bd$n2, 2] - P[bd$n1, 2]; ny <- -(P[bd$n2, 1] - P[bd$n1, 1])
  # outward normal points away from the center (0.5, 0.5)
  expect_true(all(nx * (mx - 0.5) + ny * (my - 0.5) > 0))
})

test_that("invalid connectivity is rejected", {
  P <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  expect_error(pde_mesh(P, rbind(c(1, 2, 5))), "outside")
  expect_error(pde_mesh(P, rbind(c(1, 2, 2), c(1, 3, 4))), "repeats")
  expect_error(pde_mesh(P, rbind(c(1, 2, 3), c(2, 3, 1), c(1, 3, 4))), "duplicates")
  expect_error(pde_mesh(rbind(P, c(2, 2)), rbind(c(1, 2, 3), c(1, 3, 4))), "not used")
  expect_error(pde_mesh(rbind(c(0, 0), c(1, 0), c(2, 0)), rbind(c(1, 2, 3))), "zero area")
  expect_error(pde_mesh(P, rbind(c(1, 3, 2), c(1, 3, 4))), "clockwise")
  m <- pde_mesh(P, rbind(c(1, 3, 2), c(1, 3, 4)), fix_orientation = TRUE)
  expect_equal(m$info$reoriented, 1L)
  expect_true(all(mesh_quality(m)$area > 0))
  expect_error(pde_mesh(P, rbind(c(1.5, 2, 3))), "whole-number")
})

test_that("inconsistent boundary definitions are rejected", {
  P <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  T <- rbind(c(1, 2, 3), c(1, 3, 4))
  be <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  m <- pde_mesh(P, T, be, c("a", "b", "c", "d"))
  expect_equal(sort(boundary_markers(m)), c("a", "b", "c", "d"))
  expect_error(pde_mesh(P, T, be[1:3, ], c("a", "b", "c")), "no marker")
  expect_error(pde_mesh(P, T, rbind(be, c(1, 3)), letters[1:5]), "not on the mesh boundary")
  expect_error(pde_mesh(P, T, be, c("a", "b")), "one entry")
})

test_that("mesh quality measures are exact for known triangles", {
  P <- rbind(c(0, 0), c(1, 0), c(0.5, sqrt(3) / 2))
  q <- mesh_quality(pde_mesh(P, rbind(1:3)))
  expect_equal(q$min_angle, 60)
  expect_equal(q$shape_ratio, 1)
  expect_equal(q$area, sqrt(3) / 4)
  q2 <- mesh_quality(pde_mesh_rectangle(nx = 1))
  expect_equal(q2$min_angle, c(45, 45))
  expect_equal(q2$max_angle, c(90, 90))
})

test_that("polygon meshes cover the domain and keep side markers", {
  L <- rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
  mk <- c("bottom", "right", "inner1", "inner2", "top", "left")
  m <- pde_mesh_polygon(L, h = 0.3, markers = mk)
  expect_equal(sum(mesh_quality(m)$area), 3)
  expect_lte(mesh_size(m), 0.3)
  expect_setequal(boundary_markers(m), mk)
  expect_gt(min(mesh_quality(m)$min_angle), 10)
  bn <- boundary_nodes(m, "inner1")
  expect_true(all(abs(m$nodes[bn, 2] - 1) < 1e-12))
})

test_that("set_boundary_marker relabels edges by midpoint", {
  m <- pde_mesh_rectangle(nx = 4)
  m <- set_boundary_marker(m, "inlet", function(x, y) x < 1e-12 & y > 0.5)
  expect_equal(sum(m$boundary$marker == "inlet"), 2)
})

test_that("point location and evaluation work", {
  V <- fe_space(pde_mesh_rectangle(nx = 5))
  u <- interpolate(function(x, y) 3 * x - y + 1, V)
  pts <- cbind(runif(50), runif(50))
  expect_equal(evaluate(u, pts), 3 * pts[, 1] - pts[, 2] + 1, tolerance = 1e-12)
  g <- evaluate(u, pts, gradient = TRUE)
  expect_equal(g$dx, rep(3, 50), tolerance = 1e-12)
  expect_equal(g$dy, rep(-1, 50), tolerance = 1e-12)
  expect_true(is.na(evaluate(u, 2, 2)))
  P <- evaluation_matrix(V, pts)
  expect_equal(as.numeric(P %*% u$values), 3 * pts[, 1] - pts[, 2] + 1, tolerance = 1e-12)
  expect_error(evaluation_matrix(V, 2, 2), "outside")
  # vertices and edges are found
  expect_equal(evaluate(u, c(0, 1, 0.2), c(0, 1, 0)), c(1, 3, 1.6), tolerance = 1e-12)
})
