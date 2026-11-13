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
#' @return An object of class `pde_domain`.
#' @examples
#' pde_domain("rectangle", x = c(0, 1), y = c(0, 2))
#' L <- rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
#' pde_domain("polygon", vertices = L)
#' @export
pde_domain <- function(type = c("rectangle", "polygon"),
                       x = c(0, 1), y = c(0, 1), vertices = NULL,
                       markers = "boundary") {
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
  area2 <- sum(v[, 1] * v[nxt, 2] - v[nxt, 1] * v[, 2])
  scale <- max(abs(v))^2
  if (abs(area2) <= 1e-14 * max(scale, 1)) stopf("Polygon has zero area.")
  if (polygon_self_intersects(v)) stopf("Polygon boundary intersects itself.")
  if (area2 < 0) {
    message("Polygon vertices were clockwise; reordered counterclockwise.")
    v <- v[nv:1, , drop = FALSE]
    # side i of the reversed polygon joins old vertices nv-i+1 and nv-i
    markers <- markers[c((nv - 1):1, nv)]
  }
  structure(list(type = "polygon", dim = 2L, vertices = v, markers = markers),
            class = "pde_domain")
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
  }
  cat("  markers:", paste(unique(x$markers), collapse = ", "), "\n")
  invisible(x)
}
