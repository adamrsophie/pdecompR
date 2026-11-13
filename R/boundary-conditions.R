#' Boundary conditions
#'
#' Boundary conditions are selected by boundary marker, never by coordinate
#' comparisons.
#'
#' * `dirichlet_bc(value, marker)`: `u = value` on the marked boundary.
#'   Enforced strongly by eliminating the boundary degrees of freedom and
#'   moving their contribution to the right-hand side, which keeps symmetric
#'   positive definite systems symmetric positive definite.
#' * `neumann_bc(flux, marker)`: `k du/dn = flux` (outward normal). Adds
#'   `int flux v ds` to the load.
#' * `robin_bc(alpha, value, marker)`: `k du/dn + alpha u = value`. Adds
#'   `int alpha u v ds` to the bilinear form and `int value v ds` to the load.
#'
#' `value`, `flux` and `alpha` may be numbers, `function(x, y)` or, for
#' time-dependent problems, `function(x, y, t)`.
#'
#' Dirichlet data are imposed by nodal interpolation at the boundary degrees
#' of freedom. When several Dirichlet conditions share a degree of freedom
#' (for example a corner), the condition listed last wins. A degree of
#' freedom with a Dirichlet condition ignores Neumann and Robin contributions.
#'
#' @param value Boundary value.
#' @param flux Outward flux `k du/dn`.
#' @param alpha Robin coefficient (non-negative).
#' @param marker Boundary markers; `"boundary"` selects the whole boundary.
#' @return An object of class `pde_bc`.
#' @examples
#' dirichlet_bc(0, "boundary")
#' neumann_bc(function(x, y) x, "top")
#' robin_bc(alpha = 2, value = 1, marker = "right")
#' @export
dirichlet_bc <- function(value = 0, marker = "boundary") {
  check_bc_value(value, "value")
  structure(list(type = "dirichlet", value = value, marker = as.character(marker)),
            class = "pde_bc")
}

#' @rdname dirichlet_bc
#' @export
neumann_bc <- function(flux = 0, marker = "boundary") {
  check_bc_value(flux, "flux")
  structure(list(type = "neumann", value = flux, marker = as.character(marker)),
            class = "pde_bc")
}

#' @rdname dirichlet_bc
#' @export
robin_bc <- function(alpha, value = 0, marker = "boundary") {
  check_bc_value(alpha, "alpha")
  check_bc_value(value, "value")
  if (is.numeric(alpha) && alpha < 0) stopf("`alpha` must be non-negative.")
  structure(list(type = "robin", alpha = alpha, value = value,
                 marker = as.character(marker)), class = "pde_bc")
}

check_bc_value <- function(v, name) {
  ok <- (is.numeric(v) && length(v) == 1L && is.finite(v)) || is.function(v)
  if (!ok) stopf("`%s` must be a finite number or a function.", name)
}

#' @export
print.pde_bc <- function(x, ...) {
  cat("<pde_bc>", x$type, "on", paste(x$marker, collapse = ", "), "\n")
  invisible(x)
}

normalize_bcs <- function(bcs) {
  if (is.null(bcs)) return(list())
  if (inherits(bcs, "pde_bc")) return(list(bcs))
  if (!is.list(bcs) || !all(vapply(bcs, inherits, logical(1), "pde_bc"))) {
    stopf("Boundary conditions must be `pde_bc` objects or a list of them.")
  }
  unname(bcs)
}

check_bc_markers <- function(bcs, mesh) {
  for (bc in bcs) boundary_rows(mesh, bc$marker)
  invisible(TRUE)
}

# Dirichlet degrees of freedom and values (later conditions override).
dirichlet_data <- function(space, bcs, t = NULL) {
  vals <- rep(NA_real_, space$ndof)
  for (bc in bcs) {
    if (bc$type != "dirichlet") next
    d <- boundary_dofs(space, bc$marker)
    xy <- space$dof_coords[d, , drop = FALSE]
    vals[d] <- eval_spatial(bc$value, xy[, 1], xy[, 2], t, "Dirichlet value")
  }
  dofs <- which(!is.na(vals))
  list(dofs = dofs, values = vals[dofs])
}

# Boundary matrix and load contributions of Neumann and Robin conditions.
natural_bc_terms <- function(space, bcs, t = NULL, quad_degree = NULL) {
  A <- NULL
  b <- numeric(space$ndof)
  for (bc in bcs) {
    if (bc$type == "neumann") {
      b <- b + assemble_boundary_vector(space, bc$value, bc$marker, t, quad_degree)
    } else if (bc$type == "robin") {
      R <- assemble_boundary_matrix(space, bc$alpha, bc$marker, t, quad_degree)
      A <- if (is.null(A)) R else A + R
      b <- b + assemble_boundary_vector(space, bc$value, bc$marker, t, quad_degree)
    }
  }
  list(A = A, b = b)
}
