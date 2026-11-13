#' Numerical quadrature rules
#'
#' `quadrature_interval()` returns an `n`-point Gauss-Legendre rule on the
#' reference interval `[0, 1]`, computed with the Golub-Welsch algorithm. It
#' integrates polynomials of degree `2 n - 1` exactly.
#'
#' `quadrature_triangle()` returns a rule on the reference triangle with
#' vertices `(0, 0)`, `(1, 0)` and `(0, 1)` (area 1/2) that is exact for all
#' polynomials of total degree `degree`:
#'
#' * degree 0 or 1: the centroid rule (1 point);
#' * degree 2: the three-point interior rule (Strang and Fix);
#' * degree 3 to 5: the seven-point rule of Radon (Stroud T2:5-1);
#' * higher degrees: a collapsed (Duffy) tensor product of Gauss-Legendre
#'   rules with positive weights and interior points.
#'
#' Rules are mapped to a physical triangle by the affine map
#' `x = x1 + B (xi, eta)`; weights are multiplied by `|det B|`.
#'
#' @param n Number of Gauss points.
#' @param degree Polynomial degree to integrate exactly.
#' @return A list with `points` (matrix), `weights` and `degree`, of class
#'   `fe_quadrature`.
#' @references
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
#' rules. *Mathematics of Computation*, 23, 221-230.
#'
#' Stroud, A. H. (1971). *Approximate Calculation of Multiple Integrals*.
#' Prentice-Hall.
#' @examples
#' q <- quadrature_triangle(4)
#' sum(q$weights)                               # area of the reference triangle
#' sum(q$weights * q$points[, 1]^2 * q$points[, 2]^2)  # exact: 1/180
#' @export
quadrature_interval <- function(n) {
  n <- check_count(n, "n")
  if (n == 1L) {
    t <- 0; w <- 2
  } else {
    k <- seq_len(n - 1L)
    off <- k / sqrt(4 * k^2 - 1)
    J <- matrix(0, n, n)
    J[cbind(k, k + 1L)] <- off
    J[cbind(k + 1L, k)] <- off
    ev <- eigen(J, symmetric = TRUE)
    o <- order(ev$values)
    t <- ev$values[o]
    w <- 2 * ev$vectors[1, o]^2
  }
  structure(list(points = matrix((t + 1) / 2, ncol = 1), weights = w / 2,
                 degree = 2L * n - 1L), class = "fe_quadrature")
}

#' @rdname quadrature_interval
#' @export
quadrature_triangle <- function(degree) {
  degree <- check_count(degree, "degree", min = 0L)
  if (degree <= 1L) {
    pts <- matrix(c(1, 1) / 3, 1); w <- 1 / 2; d <- 1L
  } else if (degree == 2L) {
    pts <- rbind(c(1, 1) / 6, c(2 / 3, 1 / 6), c(1 / 6, 2 / 3)); w <- rep(1 / 6, 3); d <- 2L
  } else if (degree <= 5L) {
    a <- (6 - sqrt(15)) / 21; b <- (6 + sqrt(15)) / 21
    pts <- rbind(c(1, 1) / 3, c(a, a), c(1 - 2 * a, a), c(a, 1 - 2 * a),
                 c(b, b), c(1 - 2 * b, b), c(b, 1 - 2 * b))
    w <- c(9 / 80, rep((155 - sqrt(15)) / 2400, 3), rep((155 + sqrt(15)) / 2400, 3))
    d <- 5L
  } else {
    gu <- quadrature_interval(ceiling((degree + 2) / 2))
    gv <- quadrature_interval(ceiling((degree + 1) / 2))
    u <- rep(gu$points[, 1], each = length(gv$weights))
    v <- rep(gv$points[, 1], times = length(gu$weights))
    w <- rep(gu$weights, each = length(gv$weights)) * rep(gv$weights, times = length(gu$weights)) * (1 - u)
    pts <- cbind(u, (1 - u) * v)
    d <- degree
  }
  dimnames(pts) <- NULL
  structure(list(points = pts, weights = w, degree = d), class = "fe_quadrature")
}

#' @export
print.fe_quadrature <- function(x, ...) {
  dim <- ncol(x$points)
  cat("<fe_quadrature>", if (dim == 1) "interval" else "triangle",
      "rule with", length(x$weights), "points, exact to degree", x$degree, "\n")
  invisible(x)
}

#' Integrate a function over a mesh
#'
#' Applies a triangle quadrature rule element by element.
#'
#' @param mesh A [pde_mesh()].
#' @param f A vectorized function `function(x, y)`.
#' @param degree Exactness degree of the quadrature rule.
#' @return The approximate integral.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 2)
#' integrate_mesh(mesh, function(x, y) x * y)   # exact: 1/4
#' @export
integrate_mesh <- function(mesh, f, degree = 4) {
  stopifnot(inherits(mesh, "pde_mesh"))
  q <- quadrature_triangle(degree)
  g <- affine_geometry(mesh, q)
  val <- matrix(eval_spatial(f, as.vector(g$x), as.vector(g$y), name = "f"), nrow(g$x))
  sum(val %*% q$weights * g$det)
}

# Affine maps of all triangles evaluated at reference points.
affine_geometry <- function(mesh, quad) {
  P <- mesh$nodes; T <- mesh$triangles
  x1 <- P[T[, 1], 1]; y1 <- P[T[, 1], 2]
  b11 <- P[T[, 2], 1] - x1; b12 <- P[T[, 3], 1] - x1
  b21 <- P[T[, 2], 2] - y1; b22 <- P[T[, 3], 2] - y1
  det <- b11 * b22 - b12 * b21
  xi <- quad$points[, 1]; eta <- quad$points[, 2]
  list(det = det, b11 = b11, b12 = b12, b21 = b21, b22 = b22,
       x = x1 + outer(b11, xi) + outer(b12, eta),
       y = y1 + outer(b21, xi) + outer(b22, eta))
}
