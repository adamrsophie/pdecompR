#' Computational domains
#'
#' `pde_domain()` describes the geometry on which a PDE is posed. Domains are
#' pure geometric descriptions; a mesh is built from a domain with
#' [pde_mesh_rectangle()] or [pde_mesh_polygon()]. One-dimensional and
#' three-dimensional domains are planned but not yet available.
#'
#' Supported types:
#'
#' * `"rectangle"`: the open rectangle `(x[1], x[2]) x (y[1], y[2])`. Boundary
#'   markers are `"bottom"`, `"right"`, `"top"` and `"left"`.
#' * `"polygon"`: a simple polygon given by its vertices in order. Side `i`
#'   joins vertex `i` to vertex `i + 1` (and the last vertex to the first).
#'   `markers` names each side; by default all sides are `"boundary"`.
#'
#' Every domain also answers to the marker `"boundary"`, which selects the
#' whole boundary.
#'
#' Validation rejects non-finite coordinates, zero-length sides, repeated
#' vertices, zero area and self-intersecting polygons. A clockwise polygon is
#' reversed to counterclockwise order (its side markers are reversed with it)
#' and a message reports the change.
#'
#' @param type One of `"rectangle"` or `"polygon"`.
#' @param x For rectangles, the two endpoints in the x direction.
#' @param y For rectangles, the two endpoints in the y direction.
#' @param vertices For polygons, a two-column numeric matrix of vertices.
#' @param markers For polygons, a character vector naming each side, or a
#'   single name used for every side.
#' @param holes For polygons, an optional list of vertex matrices, each a
#'   simple polygon strictly inside the outer polygon; holes must not touch
#'   each other or the outer boundary.
#' @param hole_markers Marker of the sides of each hole: one name for all
#'   holes or one per hole.
#' @return An object of class `pde_domain`.
#' @examples
#' pde_domain("rectangle", x = c(0, 1), y = c(0, 2))
#' L <- rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
#' pde_domain("polygon", vertices = L)
#' @export
pde_domain <- function(type = c("rectangle", "polygon"),
                       x = c(0, 1), y = c(0, 1), vertices = NULL,
                       markers = "boundary", holes = NULL, hole_markers = "hole") {
  type <- match.arg(type)
  if (type == "rectangle") {
    check_endpoints(x, "x")
    check_endpoints(y, "y")
    v <- rbind(c(x[1], y[1]), c(x[2], y[1]), c(x[2], y[2]), c(x[1], y[2]))
    return(structure(list(type = "rectangle", dim = 2L, x = as.numeric(x),
                          y = as.numeric(y), vertices = v,
                          markers = c("bottom", "right", "top", "left")),
                     class = "pde_domain"))
  }
  if (is.null(vertices)) stopf("A polygon domain needs `vertices`.")
  v <- as.matrix(vertices)
  storage.mode(v) <- "double"
  if (ncol(v) != 2L || nrow(v) < 3L) {
    stopf("`vertices` must be a matrix with two columns and at least three rows.")
  }
  if (any(!is.finite(v))) stopf("`vertices` contains non-finite coordinates.")
  if (nrow(v) > 3L && all(v[1, ] == v[nrow(v), ])) {
    v <- v[-nrow(v), , drop = FALSE]
  }
  nv <- nrow(v)
  if (length(markers) == 1L) markers <- rep(markers, nv)
  if (!is.character(markers) || length(markers) != nv || anyNA(markers)) {
    stopf("`markers` must be one name or one name per polygon side (%d).", nv)
  }
  if (anyDuplicated(paste(v[, 1], v[, 2]))) stopf("Polygon has repeated vertices.")
  nxt <- c(2:nv, 1L)
  len <- sqrt(rowSums((v[nxt, , drop = FALSE] - v)^2))
  if (any(len <= 0)) stopf("Polygon has a zero-length side.")
  # Shoelace formula on recentered coordinates, compared with the squared
  # extent, so the test does not depend on units or on the origin.
  vc <- sweep(v, 2, colMeans(v))
  area2 <- sum(vc[, 1] * vc[nxt, 2] - vc[nxt, 1] * vc[, 2])
  scale <- max(apply(v, 2, function(z) diff(range(z))))^2
  if (nv > 3L && polygon_self_intersects(v)) stopf("Polygon boundary intersects itself.")
  if (abs(area2) <= 1e-12 * scale) stopf("Polygon has zero area.")
  if (area2 < 0) {
    message("Polygon vertices were clockwise; reordered counterclockwise.")
    v <- v[nv:1, , drop = FALSE]
    # side i of the reversed polygon joins old vertices nv-i+1 and nv-i
    markers <- markers[c((nv - 1):1, nv)]
  }
  hl <- NULL; hm <- NULL
  if (!is.null(holes)) {
    if (is.matrix(holes)) holes <- list(holes)
    if (!is.list(holes) || !length(holes)) stopf("`holes` must be a list of vertex matrices.")
    if (length(hole_markers) == 1L) hole_markers <- rep(hole_markers, length(holes))
    if (!is.character(hole_markers) || length(hole_markers) != length(holes) || anyNA(hole_markers) ||
        any(!nzchar(hole_markers))) {
      stopf("`hole_markers` must be one name or one name per hole.")
    }
    hl <- lapply(holes, function(h) suppressMessages(pde_domain("polygon", vertices = h))$vertices)
    loops <- c(list(v), hl)
    for (i in seq_along(hl)) {
      if (!all(points_in_polygon(hl[[i]], v))) stopf("Hole %d is not strictly inside the outer polygon.", i)
      for (j in seq_along(hl)) {
        if (i != j && any(points_in_polygon(hl[[i]], hl[[j]]))) stopf("Holes %d and %d overlap.", i, j)
      }
    }
    for (i in seq_along(loops)) for (j in seq_along(loops)) {
      if (j > i && loops_intersect(loops[[i]], loops[[j]])) {
        stopf("The boundary of %s touches or crosses the boundary of hole %d.",
              if (i == 1L) "the outer polygon" else sprintf("hole %d", i - 1L), j - 1L)
      }
    }
    hm <- hole_markers
  }
  structure(list(type = "polygon", dim = 2L, vertices = v, markers = markers, holes = hl,
                 hole_markers = hm), class = "pde_domain")
}

# Ray casting; points on the boundary count as outside.
points_in_polygon <- function(p, v) {
  p <- matrix(p, ncol = 2); nv <- nrow(v); j <- c(nv, seq_len(nv - 1L))
  vapply(seq_len(nrow(p)), function(k) {
    x <- p[k, 1]; y <- p[k, 2]
    xi <- v[, 1]; yi <- v[, 2]; xj <- v[j, 1]; yj <- v[j, 2]
    cr <- (yi > y) != (yj > y)
    xs <- xj + (y - yj) * (xi - xj) / (yi - yj)
    on <- abs((xi - xj) * (y - yj) - (yi - yj) * (x - xj)) <= 1e-12 * max(1, abs(v)) &
      x >= pmin(xi, xj) & x <= pmax(xi, xj) & y >= pmin(yi, yj) & y <= pmax(yi, yj)
    !any(on) && (sum(cr & x < xs) %% 2L == 1L)
  }, logical(1))
}

# Do any sides of two closed polygons touch or cross?
loops_intersect <- function(a, b) {
  na <- nrow(a); nb <- nrow(b)
  for (i in seq_len(na)) for (j in seq_len(nb)) {
    if (segments_touch(a[i, ], a[i %% na + 1L, ], b[j, ], b[j %% nb + 1L, ])) return(TRUE)
  }
  FALSE
}

segments_touch <- function(a, b, c, d) {
  orient <- function(p, q, r) (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
  on_seg <- function(p, q, r) min(p[1], q[1]) <= r[1] && r[1] <= max(p[1], q[1]) &&
    min(p[2], q[2]) <= r[2] && r[2] <= max(p[2], q[2])
  o1 <- orient(a, b, c); o2 <- orient(a, b, d); o3 <- orient(c, d, a); o4 <- orient(c, d, b)
  if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) && ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))) return(TRUE)
  (o1 == 0 && on_seg(a, b, c)) || (o2 == 0 && on_seg(a, b, d)) ||
    (o3 == 0 && on_seg(c, d, a)) || (o4 == 0 && on_seg(c, d, b))
}

check_endpoints <- function(x, name) {
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x))) {
    stopf("`%s` must be two finite numbers.", name)
  }
  if (x[2] <= x[1]) stopf("`%s[2]` must be larger than `%s[1]`.", name, name)
}

# Proper intersection test between non-adjacent polygon sides.
polygon_self_intersects <- function(v) {
  nv <- nrow(v)
  if (nv < 4L) return(FALSE)
  orient <- function(a, b, c) {
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
  }
  on_seg <- function(a, b, c) {
    min(a[1], b[1]) <= c[1] && c[1] <= max(a[1], b[1]) &&
      min(a[2], b[2]) <= c[2] && c[2] <= max(a[2], b[2])
  }
  for (i in seq_len(nv)) {
    a <- v[i, ]; b <- v[i %% nv + 1L, ]
    for (j in seq_len(nv)) {
      if (j <= i) next
      if (j == i + 1L || (i == 1L && j == nv)) next
      c <- v[j, ]; d <- v[j %% nv + 1L, ]
      o1 <- orient(a, b, c); o2 <- orient(a, b, d)
      o3 <- orient(c, d, a); o4 <- orient(c, d, b)
      if (((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) &&
          ((o3 > 0 && o4 < 0) || (o3 < 0 && o4 > 0))) return(TRUE)
      if (o1 == 0 && on_seg(a, b, c)) return(TRUE)
      if (o2 == 0 && on_seg(a, b, d)) return(TRUE)
      if (o3 == 0 && on_seg(c, d, a)) return(TRUE)
      if (o4 == 0 && on_seg(c, d, b)) return(TRUE)
    }
  }
  FALSE
}

#' @export
print.pde_domain <- function(x, ...) {
  cat("<pde_domain>", x$type, "\n")
  if (x$type == "rectangle") {
    cat("  (", x$x[1], ", ", x$x[2], ") x (", x$y[1], ", ", x$y[2], ")\n", sep = "")
  } else {
    cat("  vertices:", nrow(x$vertices), "\n")
    if (length(x$holes)) cat("  holes:", length(x$holes), "\n")
  }
  cat("  markers:", paste(unique(x$markers), collapse = ", "), "\n")
  invisible(x)
}
