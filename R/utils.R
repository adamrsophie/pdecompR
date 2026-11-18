#' @keywords internal
#' @importFrom Matrix sparseMatrix t crossprod tcrossprod Diagonal Cholesky solve diag drop0
#' @importFrom methods as is
#' @importFrom stats optim setNames
#' @importFrom graphics plot polygon segments lines points legend par contour image axis box text abline rect
#' @importFrom grDevices hcl.colors
#' @importFrom utils head tail
"_PACKAGE"

# Internal helpers ---------------------------------------------------------

stopf <- function(...) stop(sprintf(...), call. = FALSE)

warnf <- function(...) warning(sprintf(...), call. = FALSE)

# Error for non-finite values produced by a user callback. Newton line
# searches catch this class (and only this class) at trial points and take a
# shorter step; any other error still propagates.
stop_nonfinite <- function(...) {
  stop(structure(class = c("pde_nonfinite", "error", "condition"),
                 list(message = sprintf(...), call = NULL)))
}

# Evaluate f() at a Newton trial point: NULL when it produced non-finite
# values. Warnings from rejected trial points (for example log of a negative
# number) are muffled; the accepted state is evaluated again without this.
trial_eval <- function(f) {
  tryCatch(suppressWarnings(f()), pde_nonfinite = function(e) NULL)
}

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

# Shared validators, one per concept.

# Point coordinates: equal-length finite numeric vectors (or a two-column
# matrix in `x` with `y = NULL`). Returns list(x, y).
check_points <- function(x, y = NULL, name = "points") {
  if (is.null(y)) {
    x <- as.matrix(x)
    if (ncol(x) == 1L && nrow(x) == 2L) x <- t(x)
    if (ncol(x) != 2L) stopf("Give %s as `x` and `y` vectors or as a two-column matrix.", name)
    y <- x[, 2]; x <- x[, 1]
  }
  if (!is.numeric(x) || !is.numeric(y)) stopf("The coordinates of %s must be numeric.", name)
  if (length(x) != length(y)) stopf("`x` and `y` must have the same length.")
  if (any(!is.finite(x)) || any(!is.finite(y))) stopf("The coordinates of %s must be finite.", name)
  list(x = as.numeric(x), y = as.numeric(y))
}

# Element (or other) indices: whole numbers in 1..n, or a logical vector of
# length n. Returns unique integer indices.
check_indices <- function(i, n, name) {
  if (is.logical(i)) {
    if (length(i) != n || anyNA(i)) stopf("`%s` must have one TRUE/FALSE value per element.", name)
    return(which(i))
  }
  if (!is.numeric(i) || anyNA(i) || any(i != round(i)) || any(i < 1 | i > n)) {
    stopf("`%s` must contain whole-number indices between 1 and %d.", name, n)
  }
  unique(as.integer(i))
}

# Boundary marker or region name.
check_name <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stopf("`%s` must be one nonempty name.", name)
  }
  x
}

# Logical selection returned by a user predicate.
check_selection <- function(sel, n, what) {
  if (!is.logical(sel) || length(sel) != n || anyNA(sel)) {
    stopf("`which` must return one TRUE/FALSE value (no NA) per %s.", what)
  }
  sel
}

# Positive finite tolerance.
check_tol <- function(x, name) {
  check_number(x, name, positive = TRUE)
  x
}

# Evaluate a user supplied function of space (and optionally time) at points.
# `f` may be a number, a vector of nodal length handled elsewhere, or a
# function(x, y) / function(x, y, t). Always returns a numeric vector of
# length nrow(points).
eval_spatial <- function(f, x, y, t = NULL, name = "function") {
  n <- length(x)
  if (is.null(f)) return(numeric(n))
  if (is.numeric(f)) {
    if (length(f) == 1L && is.finite(f)) return(rep(as.numeric(f), n))
    stopf("`%s` must be a single number or a function of (x, y).", name)
  }
  if (inherits(f, "fe_function")) {
    val <- evaluate(f, x, y)
    if (anyNA(val)) stopf("`%s` (a finite-element function) does not cover all points.", name)
    return(val)
  }
  if (!is.function(f)) stopf("`%s` must be a number or a function.", name)
  td <- is_time_function(f)
  if (td && is.null(t)) {
    stopf("`%s` depends on time (third argument `t`) but the problem is stationary.", name)
  }
  val <- if (td) f(x, y, t) else f(x, y)
  if (!is.numeric(val)) stopf("`%s` must return numeric values.", name)
  if (length(val) == 1L) val <- rep(val, n)
  if (length(val) != n) {
    stopf("`%s` returned %d values for %d points; it must be vectorized.",
          name, length(val), n)
  }
  if (!all(is.finite(val))) stopf("`%s` returned non-finite values (NA, NaN or Inf).", name)
  as.numeric(val)
}

# A function is time dependent when its third formal argument is named `t`;
# extra arguments with other names (for example `amp = 5`) keep their defaults.
is_time_function <- function(f) {
  # A tensor depends on time when any of its entries does.
  if (inherits(f, "pde_tensor")) return(any(vapply(f[c("xx", "xy", "yy")], is_time_function, logical(1))))
  is.function(f) && length(formals(f)) >= 3L && identical(names(formals(f))[3], "t")
}

fmt_num <- function(x, digits = 4) formatC(as.numeric(x), digits = digits, format = "g")

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
