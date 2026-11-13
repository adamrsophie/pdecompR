test_that("rectangle domains validate endpoints", {
  d <- pde_domain("rectangle", x = c(0, 2), y = c(-1, 1))
  expect_s3_class(d, "pde_domain")
  expect_equal(d$markers, c("bottom", "right", "top", "left"))
  expect_error(pde_domain("rectangle", x = c(1, 0)), "larger")
  expect_error(pde_domain("rectangle", x = c(0, NA)), "finite")
})

test_that("polygon validation catches degenerate input", {
  sq <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  expect_s3_class(pde_domain("polygon", vertices = sq), "pde_domain")
  expect_error(pde_domain("polygon", vertices = rbind(c(0, 0), c(1, 0), c(2, 0))), "zero area")
  expect_error(pde_domain("polygon", vertices = rbind(c(0, 0), c(1, 0), c(1, 0), c(0, 1))), "repeated")
  bow <- rbind(c(0, 0), c(1, 1), c(1, 0), c(0, 1))
  expect_error(pde_domain("polygon", vertices = bow), "intersects")
  expect_error(pde_domain("polygon", vertices = rbind(c(0, 0), c(Inf, 0), c(0, 1))), "non-finite")
  expect_error(pde_domain("polygon", vertices = sq, markers = c("a", "b")), "one name")
})

test_that("clockwise polygons are reoriented with their markers", {
  cw <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0))
  expect_message(d <- pde_domain("polygon", vertices = cw, markers = c("left", "top", "right", "bottom")),
                 "counterclockwise")
  v <- d$vertices
  nxt <- c(2:4, 1)
  area2 <- sum(v[, 1] * v[nxt, 2] - v[nxt, 1] * v[, 2])
  expect_gt(area2, 0)
  # each side keeps its marker: check the side on y = 0 is "bottom"
  mid <- (v + v[nxt, ]) / 2
  expect_equal(d$markers[which(abs(mid[, 2]) < 1e-12)], "bottom")
  expect_equal(d$markers[which(abs(mid[, 1]) < 1e-12)], "left")
})
