#' Lagrange finite elements on triangles
#'
#' `fe_element()` describes the continuous Lagrange element of degree 1 (P1)
#' or 2 (P2) on the reference triangle `(0, 0)`, `(1, 0)`, `(0, 1)`. With
#' barycentric coordinates `l1 = 1 - xi - eta`, `l2 = xi`, `l3 = eta`:
#'
#' * P1: `phi_i = l_i`, `i = 1, 2, 3`, with nodes at the vertices.
#' * P2: `phi_i = l_i (2 l_i - 1)` at the vertices and
#'   `phi_4 = 4 l1 l2`, `phi_5 = 4 l2 l3`, `phi_6 = 4 l3 l1` at the midpoints
#'   of local edges 1, 2 and 3.
#'
#' Both satisfy the nodal property `phi_i(x_j) = delta_ij`.
#'
#' @param degree Polynomial degree, 1 or 2.
#' @return An object of class `fe_element` with fields `degree`, `ndof`,
#'   `nodes` (reference coordinates of the local nodes) and functions
#'   `basis(xi, eta)` and `grad(xi, eta)`. `basis` returns a matrix with one
#'   row per point and one column per local basis function; `grad` returns a
#'   list with the `xi` and `eta` derivatives in the same layout.
#' @examples
#' el <- fe_element(1)
#' el$basis(el$nodes[, 1], el$nodes[, 2])   # identity matrix
#' @export
fe_element <- function(degree = 1) {
  degree <- check_count(degree, "degree")
  if (degree == 1L) {
    basis <- function(xi, eta) cbind(1 - xi - eta, xi, eta)
    grad <- function(xi, eta) {
      n <- length(xi)
      list(xi = cbind(rep(-1, n), rep(1, n), rep(0, n)),
           eta = cbind(rep(-1, n), rep(0, n), rep(1, n)))
    }
    nodes <- rbind(c(0, 0), c(1, 0), c(0, 1))
  } else if (degree == 2L) {
    basis <- function(xi, eta) {
      l1 <- 1 - xi - eta; l2 <- xi; l3 <- eta
      cbind(l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1),
            4 * l1 * l2, 4 * l2 * l3, 4 * l3 * l1)
    }
    grad <- function(xi, eta) {
      l1 <- 1 - xi - eta; l2 <- xi; l3 <- eta
      z <- rep(0, length(xi))
      list(xi = cbind(-(4 * l1 - 1), 4 * l2 - 1, z,
                      4 * (l1 - l2), 4 * l3, -4 * l3),
           eta = cbind(-(4 * l1 - 1), z, 4 * l3 - 1,
                       -4 * l2, 4 * l2, 4 * (l1 - l3)))
    }
    nodes <- rbind(c(0, 0), c(1, 0), c(0, 1), c(0.5, 0), c(0.5, 0.5), c(0, 0.5))
  } else {
    stopf("Only Lagrange elements of degree 1 and 2 are implemented.")
  }
  structure(list(family = "Lagrange", degree = degree, ndof = nrow(nodes),
                 nodes = nodes, basis = basis, grad = grad),
            class = "fe_element")
}

#' @export
print.fe_element <- function(x, ...) {
  cat("<fe_element> continuous Lagrange P", x$degree, " on triangles, ",
      x$ndof, " local basis functions\n", sep = "")
  invisible(x)
}

#' Local P1 element matrices
#'
#' Closed-form P1 stiffness and mass matrices of a single triangle, used as
#' independent references for the quadrature-based assembly.
#'
#' The stiffness matrix is `K_e = |T| G G^T` where the rows of `G` are the
#' constant gradients of the three basis functions; the mass matrix is
#' `M_e = |T| / 12 * (1 + I)` (2 on the diagonal, 1 elsewhere).
#'
#' @param vertices A 3 by 2 matrix of vertex coordinates (counterclockwise).
#' @param k Constant diffusion coefficient.
#' @return A 3 by 3 matrix.
#' @examples
#' v <- rbind(c(0, 0), c(1, 0), c(0, 1))
#' element_stiffness_p1(v)
#' element_mass_p1(v)
#' @export
element_stiffness_p1 <- function(vertices, k = 1) {
  v <- as.matrix(vertices)
  B <- cbind(v[2, ] - v[1, ], v[3, ] - v[1, ])
  area <- det(B) / 2
  if (area <= 0) stopf("Vertices must be counterclockwise and non-degenerate.")
  G <- rbind(c(-1, -1), c(1, 0), c(0, 1)) %*% solve(B)
  k * area * G %*% t(G)
}

#' @rdname element_stiffness_p1
#' @export
element_mass_p1 <- function(vertices) {
  v <- as.matrix(vertices)
  area <- det(cbind(v[2, ] - v[1, ], v[3, ] - v[1, ])) / 2
  if (area <= 0) stopf("Vertices must be counterclockwise and non-degenerate.")
  area / 12 * (matrix(1, 3, 3) + diag(3))
}
