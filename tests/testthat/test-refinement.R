test_that("uniform refinement quarters elements and keeps markers", {
  m <- pde_mesh_rectangle(nx = 2, ny = 3)
  r <- refine_uniform(m)
  expect_equal(nrow(r$triangles), 4 * nrow(m$triangles))
  expect_equal(nrow(r$boundary), 2 * nrow(m$boundary))
  expect_equal(sum(mesh_quality(r)$area), 1)
  expect_equal(table(r$boundary$marker)[["left"]], 2 * table(m$boundary$marker)[["left"]])
  expect_equal(min(mesh_quality(r)$min_angle), min(mesh_quality(m)$min_angle))
  expect_equal(mesh_size(r), mesh_size(m) / 2)
})

test_that("newest vertex bisection is conforming and shape regular", {
  m <- pde_mesh_rectangle(nx = 4)
  for (k in 1:6) {
    # refine around the corner (0, 0)
    P <- m$nodes; T <- m$triangles
    cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
    cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
    m <- refine_mesh(m, which(cx^2 + cy^2 < 0.1))
    expect_equal(sum(mesh_quality(m)$area), 1)
    # conformity: every interior edge has two elements, Euler characteristic 1
    expect_equal(nrow(m$nodes) - nrow(m$edges) + nrow(m$triangles), 1)
  }
  expect_gt(min(mesh_quality(m)$min_angle), 20)
  expect_setequal(boundary_markers(m), c("bottom", "right", "top", "left"))
  expect_lt(min(mesh_quality(m)$h_min), 0.25 / 4)
})

test_that("refinement records parents and transfers P1 functions exactly", {
  m <- pde_mesh_rectangle(nx = 3)
  V <- fe_space(m)
  u <- interpolate(function(x, y) x^2 + y, V)
  r <- refine_mesh(m, c(1, 5))
  par <- r$info$new_nodes
  expect_equal(nrow(r$nodes), nrow(m$nodes) + nrow(par))
  ut <- transfer_function(u, r)
  new <- nrow(m$nodes) + seq_len(nrow(par))
  expect_equal(ut$values[new], (u$values[par[, 1]] + u$values[par[, 2]]) / 2, tolerance = 1e-12)
})

test_that("refinement records parent elements and inherits regions through them", {
  m <- set_region(pde_mesh_rectangle(nx = 4), "a", function(x, y) x < 0.5 & y > 0.25)
  m <- set_region(m, "b", function(x, y) x > 0.6)
  f <- refine_uniform(m); g <- refine_mesh(f, c(1, 5, 9, 30)); g2 <- refine_mesh(g, 1:20)
  for (k in list(list(f, m), list(g, f), list(g2, g))) {
    ch <- k[[1]]; pa <- k[[2]]
    P <- ch$nodes; T <- ch$triangles
    cx <- rowMeans(matrix(P[T, 1], ncol = 3)); cy <- rowMeans(matrix(P[T, 2], ncol = 3))
    expect_equal(locate_points(pa, cx, cy)$element, ch$info$parent_elements)
    expect_equal(mesh_regions(ch), mesh_regions(pa)[ch$info$parent_elements])
  }
})
