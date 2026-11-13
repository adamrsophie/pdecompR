#' @keywords internal
#' @importFrom Matrix sparseMatrix t crossprod tcrossprod Diagonal Cholesky solve diag drop0
#' @importFrom methods as is
#' @importFrom stats optim setNames
#' @importFrom graphics plot polygon segments lines points legend par contour image axis box text abline
#' @importFrom grDevices hcl.colors
#' @importFrom utils head tail
"_PACKAGE"

# Internal helpers ---------------------------------------------------------

stopf <- function(...) stop(sprintf(...), call. = FALSE)

warnf <- function(...) warning(sprintf(...), call. = FALSE)

check_number <- function(x, name, positive = FALSE, finite = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stopf("`%s` must be a single number.", name)
  }
  if (finite && !is.finite(x)) stopf("`%s` must be finite.", name)
  if (positive && x <= 0) stopf("`%s` must be positive.", name)
  invisible(x)
}

check_count <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != round(x) || x < min) {
    stopf("`%s` must be a whole number greater than or equal to %d.", name, min)
  }
  as.integer(x)
}

# Evaluate a user supplied function of space (and optionally time) at points.
# `f` may be a number, a vector of nodal length handled elsewhere, or a
# function(x, y) / function(x, y, t). Always returns a numeric vector of
# length nrow(points).
eval_spatial <- function(f, x, y, t = NULL, name = "function") {
  n <- length(x)
  if (is.null(f)) return(numeric(n))
  if (is.numeric(f)) {
    if (length(f) == 1L) return(rep(as.numeric(f), n))
    stopf("`%s` must be a single number or a function of (x, y).", name)
  }
  if (!is.function(f)) stopf("`%s` must be a number or a function.", name)
  nf <- length(formals(f))
  has_dots <- "..." %in% names(formals(f))
  val <- if (!is.null(t) && (nf >= 3L || has_dots)) f(x, y, t) else f(x, y)
  if (!is.numeric(val)) stopf("`%s` must return numeric values.", name)
  if (length(val) == 1L) val <- rep(val, n)
  if (length(val) != n) {
    stopf("`%s` returned %d values for %d points; it must be vectorized.",
          name, length(val), n)
  }
  as.numeric(val)
}

is_time_function <- function(f) {
  is.function(f) && (length(formals(f)) >= 3L || "..." %in% names(formals(f)))
}

fmt_num <- function(x, digits = 4) formatC(x, digits = digits, format = "g")

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, time = proc.time()[["elapsed"]] - t0)
}

# Signed doubled areas of triangles given coordinates.
tri_det <- function(nodes, tris) {
  x1 <- nodes[tris[, 1], 1]; y1 <- nodes[tris[, 1], 2]
  x2 <- nodes[tris[, 2], 1]; y2 <- nodes[tris[, 2], 2]
  x3 <- nodes[tris[, 3], 1]; y3 <- nodes[tris[, 3], 2]
  (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
}

edge_key <- function(a, b, n) {
  lo <- pmin(a, b); hi <- pmax(a, b)
  (as.numeric(lo) - 1) * n + as.numeric(hi)
}
