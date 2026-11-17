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

#' Periodic boundary conditions
#'
#' Identifies the degrees of freedom on boundary `to` with those on boundary
#' `from` translated by `offset`, so that `u(x + offset) = u(x)` (and, in the
#' weak form, the fluxes through the two boundaries cancel). Pairs are matched
#' by coordinates; every degree of freedom on `to` must have a partner on
#' `from`, so the two boundaries need matching nodes, as in
#' [pde_mesh_rectangle()]. Several periodic conditions (for example in `x`
#' and `y`) may be combined, and corners are identified transitively.
#' Supported for [elliptic_problem()]; with no Dirichlet condition and no
#' reaction the solution is fixed by a zero-mean gauge.
#'
#' @param from,to Boundary markers.
#' @param offset Translation `c(dx, dy)` taking `from` onto `to`.
#' @return An object of class `pde_bc`.
#' @examples
#' m <- pde_mesh_rectangle(nx = 16)
#' p <- elliptic_problem(m, source = function(x, y) sin(2 * pi * x),
#'                       bcs = list(periodic_bc("left", "right", c(1, 0)),
#'                                  dirichlet_bc(0, c("bottom", "top"))))
#' range(solve_pde(p)$u$values)
#' @export
periodic_bc <- function(from, to, offset) {
  check_name(from, "from"); check_name(to, "to")
  if (!is.numeric(offset) || length(offset) != 2L || any(!is.finite(offset)) || all(offset == 0)) {
    stopf("`offset` must be a nonzero finite translation c(dx, dy).")
  }
  structure(list(type = "periodic", from = from, to = to, offset = as.numeric(offset),
                 marker = c(from, to)), class = "pde_bc")
}

# Constraint matrix of periodic conditions: P (ndof by nr) with u = P u_r,
# and the reduced index of every dof. NULL without periodic conditions.
periodic_constraints <- function(space, bcs) {
  pb <- Filter(function(b) b$type == "periodic", bcs)
  if (!length(pb)) return(NULL)
  n <- space$ndof
  rep_of <- seq_len(n)
  X <- space$dof_coords
  scale <- max(abs(X), 1)
  for (b in pb) {
    df <- boundary_dofs(space, b$from); dt <- boundary_dofs(space, b$to)
    key <- function(xy) paste(round(xy[, 1] / scale, 9), round(xy[, 2] / scale, 9))
    kf <- key(sweep(X[df, , drop = FALSE], 2, b$offset, `+`))
    m <- match(key(X[dt, , drop = FALSE]), kf)
    if (anyNA(m)) {
      stopf("Periodic condition %s -> %s: %d degree(s) of freedom on `%s` have no partner on `%s` at offset (%s, %s).",
            b$from, b$to, sum(is.na(m)), b$to, b$from, fmt_num(b$offset[1]), fmt_num(b$offset[2]))
    }
    rep_of[dt] <- df[m]
  }
  # Follow chains (corners) to the root.
  for (k in seq_len(n)) { nx <- rep_of[rep_of]; if (identical(nx, rep_of)) break; rep_of <- nx }
  if (any(rep_of[rep_of] != rep_of)) stopf("Periodic conditions form a cycle; check the offsets.")
  roots <- sort(unique(rep_of))
  col <- match(rep_of, roots)
  list(P = Matrix::sparseMatrix(i = seq_len(n), j = col, x = 1, dims = c(n, length(roots))), col = col)
}

no_periodic <- function(bcs, what) {
  if (any(vapply(bcs, function(b) identical(b$type, "periodic"), logical(1)))) {
    stopf("Periodic conditions are supported for elliptic_problem() only, not for %s.", what)
  }
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
# With `matrix = FALSE` only the load is assembled (time stepping, where the
# Robin matrix is assembled once).
natural_bc_terms <- function(space, bcs, t = NULL, quad_degree = NULL, matrix = TRUE) {
  A <- NULL
  b <- numeric(space$ndof)
  for (bc in bcs) {
    if (!bc$type %in% c("neumann", "robin")) next
    if (matrix && bc$type == "robin") {
      R <- assemble_boundary_matrix(space, bc$alpha, bc$marker, t, quad_degree)
      A <- if (is.null(A)) R else A + R
    }
    b <- b + assemble_boundary_vector(space, bc$value, bc$marker, t, quad_degree)
  }
  list(A = A, b = b)
}
