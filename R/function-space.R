#' Finite-element function spaces
#'
#' `fe_space()` builds a continuous Lagrange space on a triangular mesh.
#'
#' @section Degrees of freedom:
#' For P1 the degrees of freedom are the mesh nodes, in node order. For P2
#' they are the nodes followed by the edge midpoints: the midpoint of mesh
#' edge `k` has index `nrow(mesh$nodes) + k`. `cell_dofs` maps each element
#' to its global degrees of freedom in local basis order.
#'
#' @param mesh A [pde_mesh()].
#' @param degree Polynomial degree, 1 or 2.
#' @param family Element family. Only `"Lagrange"` is implemented.
#' @return An object of class `fe_space` with fields `mesh`, `element`,
#'   `ndof`, `cell_dofs` and `dof_coords`.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 4), degree = 2)
#' V
#' @export
fe_space <- function(mesh, degree = 1, family = "Lagrange") {
  stopifnot(inherits(mesh, "pde_mesh"))
  if (!identical(family, "Lagrange")) {
    stopf("Only the continuous Lagrange family is implemented.")
  }
  el <- fe_element(degree)
  n <- nrow(mesh$nodes)
  if (el$degree == 1L) {
    cell_dofs <- mesh$triangles
    coords <- mesh$nodes
  } else {
    cell_dofs <- cbind(mesh$triangles, n + mesh$element_edges)
    E <- mesh$edges
    coords <- rbind(mesh$nodes,
                    (mesh$nodes[E[, 1], , drop = FALSE] + mesh$nodes[E[, 2], , drop = FALSE]) / 2)
  }
  dimnames(cell_dofs) <- NULL
  structure(list(mesh = mesh, element = el, degree = el$degree,
                 ndof = nrow(coords), cell_dofs = cell_dofs, dof_coords = coords),
            class = "fe_space")
}

#' @export
print.fe_space <- function(x, ...) {
  cat("<fe_space> continuous Lagrange P", x$degree, ", ", x$ndof,
      " degrees of freedom on ", nrow(x$mesh$triangles), " triangles\n", sep = "")
  invisible(x)
}

#' Boundary degrees of freedom
#'
#' @param space An [fe_space()].
#' @param marker Boundary markers; `"boundary"` selects the whole boundary.
#' @return Sorted integer indices of the degrees of freedom on the selected
#'   boundary edges.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 2))
#' boundary_dofs(V, "bottom")
#' @export
boundary_dofs <- function(space, marker = "boundary") {
  mesh <- space$mesh
  rows <- boundary_rows(mesh, marker)
  d <- c(mesh$boundary$n1[rows], mesh$boundary$n2[rows])
  if (space$degree == 2L) d <- c(d, nrow(mesh$nodes) + mesh$boundary$edge[rows])
  sort(unique(d))
}

#' Finite-element functions
#'
#' `fe_function()` wraps a coefficient vector in a space. `interpolate()`
#' returns the nodal (Lagrange) interpolant of a function.
#'
#' @param space An [fe_space()].
#' @param values Numeric vector of length `space$ndof`.
#' @param name Optional label used in printing and plots.
#' @param f A number or vectorized function `function(x, y)`.
#' @return An object of class `fe_function`.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 4))
#' u <- interpolate(function(x, y) x + 2 * y, V)
#' evaluate(u, 0.3, 0.4)
#' @export
fe_function <- function(space, values, name = "u") {
  stopifnot(inherits(space, "fe_space"))
  if (!is.numeric(values)) stopf("`values` must be numeric.")
  values <- as.numeric(values)
  if (length(values) != space$ndof) {
    stopf("`values` has length %d; the space has %d degrees of freedom.",
          length(values), space$ndof)
  }
  structure(list(space = space, values = values, name = name), class = "fe_function")
}

#' @rdname fe_function
#' @export
interpolate <- function(f, space, name = "u") {
  stopifnot(inherits(space, "fe_space"))
  if (inherits(f, "fe_function")) return(transfer_function(f, space$mesh, space$degree))
  xy <- space$dof_coords
  fe_function(space, eval_spatial(f, xy[, 1], xy[, 2], name = "f"), name = name)
}

#' @export
print.fe_function <- function(x, ...) {
  cat("<fe_function> ", x$name, " in P", x$space$degree, " (", x$space$ndof,
      " dofs), range [", fmt_num(min(x$values)), ", ", fmt_num(max(x$values)), "]\n",
      sep = "")
  invisible(x)
}

#' Evaluate finite-element functions at points
#'
#' `evaluate()` returns values (and optionally gradients) of a finite-element
#' function at arbitrary points. `evaluation_matrix()` returns the sparse
#' matrix `P` with `P %*% values` equal to the point values, which is the
#' observation operator used for inverse problems.
#'
#' Points are located with a uniform bucket grid. Points outside the mesh
#' give `NA` (in `evaluate()`) or an error (in `evaluation_matrix()`).
#'
#' @param u An [fe_function()].
#' @param x,y Coordinates (vectors of equal length), or a two-column matrix in
#'   `x`.
#' @param gradient Also return the gradient.
#' @param space An [fe_space()].
#' @return `evaluate()` returns a numeric vector, or with `gradient = TRUE` a
#'   data frame with columns `value`, `dx`, `dy`. `evaluation_matrix()`
#'   returns a sparse `dgCMatrix`.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 4))
#' u <- interpolate(function(x, y) x^2, V)
#' evaluate(u, c(0.5, 0.25), c(0.5, 0.5), gradient = TRUE)
#' @export
evaluate <- function(u, x, y = NULL, gradient = FALSE) {
  stopifnot(inherits(u, "fe_function"))
  pts <- check_points(x, y); x <- pts$x; y <- pts$y
  loc <- locate_points(u$space$mesh, x, y)
  ok <- !is.na(loc$element)
  val <- rep(NA_real_, length(x))
  if (gradient) { gx <- val; gy <- val }
  if (any(ok)) {
    e <- loc$element[ok]
    el <- u$space$element
    B <- el$basis(loc$xi[ok], loc$eta[ok])
    D <- u$space$cell_dofs[e, , drop = FALSE]
    U <- matrix(u$values[D], ncol = el$ndof)
    val[ok] <- rowSums(B * U)
    if (gradient) {
      G <- el$grad(loc$xi[ok], loc$eta[ok])
      gr <- physical_gradients(u$space$mesh, e, G)
      gx[ok] <- rowSums(gr$x * U); gy[ok] <- rowSums(gr$y * U)
    }
  }
  if (gradient) return(data.frame(value = val, dx = gx, dy = gy))
  val
}

#' @rdname evaluate
#' @export
evaluation_matrix <- function(space, x, y = NULL) {
  stopifnot(inherits(space, "fe_space"))
  pts <- check_points(x, y); x <- pts$x; y <- pts$y
  loc <- locate_points(space$mesh, x, y)
  if (anyNA(loc$element)) {
    stopf("%d point(s) lie outside the mesh (first: %d).",
          sum(is.na(loc$element)), which(is.na(loc$element))[1])
  }
  B <- space$element$basis(loc$xi, loc$eta)
  D <- space$cell_dofs[loc$element, , drop = FALSE]
  Matrix::sparseMatrix(i = rep(seq_along(x), ncol(B)), j = as.vector(D),
                       x = as.vector(B), dims = c(length(x), space$ndof))
}

# Physical gradients of local basis functions in elements `e`, given the
# reference gradients (list of xi and eta matrices, one row per element).
physical_gradients <- function(mesh, e, G) {
  P <- mesh$nodes; T <- mesh$triangles[e, , drop = FALSE]
  x1 <- P[T[, 1], 1]; y1 <- P[T[, 1], 2]
  b11 <- P[T[, 2], 1] - x1; b12 <- P[T[, 3], 1] - x1
  b21 <- P[T[, 2], 2] - y1; b22 <- P[T[, 3], 2] - y1
  det <- b11 * b22 - b12 * b21
  list(x = (b22 * G$xi - b21 * G$eta) / det,
       y = (-b12 * G$xi + b11 * G$eta) / det)
}

# Locate points in the mesh. Returns the element and reference coordinates.
# ponytail: uniform bucket grid; a k-d tree would help strongly graded meshes.
build_search_index <- function(mesh) {
  P <- mesh$nodes; T <- mesh$triangles
  m <- nrow(T)
  xmin <- min(P[, 1]); xmax <- max(P[, 1]); ymin <- min(P[, 2]); ymax <- max(P[, 2])
  nb <- max(1L, ceiling(sqrt(m / 2)))
  wx <- (xmax - xmin) / nb; wy <- (ymax - ymin) / nb
  bin <- function(v, lo, w) pmin(nb - 1L, pmax(0L, floor((v - lo) / w)))
  ex <- matrix(P[T, 1], m); ey <- matrix(P[T, 2], m)
  slack <- 1e-9 * max(xmax - xmin, ymax - ymin)
  bx0 <- bin(pmin(ex[, 1], ex[, 2], ex[, 3]) - slack, xmin, wx)
  bx1 <- bin(pmax(ex[, 1], ex[, 2], ex[, 3]) + slack, xmin, wx)
  by0 <- bin(pmin(ey[, 1], ey[, 2], ey[, 3]) - slack, ymin, wy)
  by1 <- bin(pmax(ey[, 1], ey[, 2], ey[, 3]) + slack, ymin, wy)
  nxk <- bx1 - bx0 + 1L; nyk <- by1 - by0 + 1L
  cnt <- nxk * nyk
  el <- rep(seq_len(m), cnt)
  off <- sequence(cnt) - 1L
  bxx <- rep(bx0, cnt) + off %% rep(nxk, cnt)
  byy <- rep(by0, cnt) + off %/% rep(nxk, cnt)
  bid <- bxx + byy * nb + 1L
  o <- order(bid)
  el <- el[o]; bid <- bid[o]
  start <- match(seq_len(nb * nb), bid)
  num <- tabulate(bid, nbins = nb * nb)
  list(nodes = P, triangles = T, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, nb = nb,
       wx = wx, wy = wy, slack = slack, el = el, start = start, num = num)
}

# Bucket index of the most recently searched meshes, keyed by their nodes
# and triangles only (markers and regions do not affect point location).
search_cache <- new.env(parent = emptyenv())
search_index <- function(mesh) {
  for (ix in search_cache$items) {
    if (identical(ix$triangles, mesh$triangles) && identical(ix$nodes, mesh$nodes)) return(ix)
  }
  ix <- build_search_index(mesh)
  search_cache$items <- c(list(ix), utils::head(search_cache$items, 3L))
  ix
}

locate_points <- function(mesh, x, y, tol = 1e-10) {
  ix <- search_index(mesh)
  P <- ix$nodes; T <- ix$triangles; np <- length(x)
  xmin <- ix$xmin; xmax <- ix$xmax; ymin <- ix$ymin; ymax <- ix$ymax
  nb <- ix$nb; wx <- ix$wx; wy <- ix$wy; slack <- ix$slack
  el <- ix$el; start <- ix$start; num <- ix$num
  bin <- function(v, lo, w) pmin(nb - 1L, pmax(0L, floor((v - lo) / w)))
  inside <- x >= xmin - slack & x <= xmax + slack & y >= ymin - slack & y <= ymax + slack
  pb <- rep(NA_integer_, np)
  pb[inside] <- bin(x[inside], xmin, wx) + bin(y[inside], ymin, wy) * nb + 1L
  res_e <- rep(NA_integer_, np); res_xi <- rep(NA_real_, np); res_eta <- rep(NA_real_, np)
  cand_n <- ifelse(is.na(pb), 0L, num[pb])
  if (sum(cand_n) == 0L) return(list(element = res_e, xi = res_xi, eta = res_eta))
  pid <- rep(seq_len(np), cand_n)
  cid <- el[rep(start[pb[cand_n > 0]], cand_n[cand_n > 0]) + sequence(cand_n[cand_n > 0]) - 1L]
  x1 <- P[T[cid, 1], 1]; y1 <- P[T[cid, 1], 2]
  b11 <- P[T[cid, 2], 1] - x1; b12 <- P[T[cid, 3], 1] - x1
  b21 <- P[T[cid, 2], 2] - y1; b22 <- P[T[cid, 3], 2] - y1
  det <- b11 * b22 - b12 * b21
  dx <- x[pid] - x1; dy <- y[pid] - y1
  xi <- (b22 * dx - b12 * dy) / det
  eta <- (-b21 * dx + b11 * dy) / det
  ok <- xi >= -tol & eta >= -tol & xi + eta <= 1 + tol
  hit <- which(ok)
  hit <- hit[!duplicated(pid[hit])]
  res_e[pid[hit]] <- cid[hit]
  res_xi[pid[hit]] <- xi[hit]; res_eta[pid[hit]] <- eta[hit]
  list(element = res_e, xi = res_xi, eta = res_eta)
}

#' Transfer a finite-element function to another mesh
#'
#' Interpolates `u` into the Lagrange space of the given degree on `mesh` by
#' evaluating `u` at the new degrees of freedom. When `mesh` was produced from
#' `u`'s mesh by [refine_uniform()] or [refine_mesh()] and `u` is P1, the
#' transfer is exact (new nodes are edge midpoints).
#'
#' @param u An [fe_function()].
#' @param mesh Target [pde_mesh()].
#' @param degree Degree of the target space.
#' @return An [fe_function()] on `mesh`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 2)
#' u <- interpolate(function(x, y) x + y, fe_space(mesh))
#' v <- transfer_function(u, refine_uniform(mesh))
#' evaluate(v, 0.3, 0.4)
#' @export
transfer_function <- function(u, mesh, degree = u$space$degree) {
  stopifnot(inherits(u, "fe_function"))
  V <- fe_space(mesh, degree)
  # One refinement step of u's mesh: P1 values at new nodes are the averages
  # of the two parent edge endpoints (exact, and no point location).
  E <- if (degree == 1L && u$space$degree == 1L) refinement_edges(u$space$mesh, mesh) else NULL
  if (!is.null(E)) return(fe_function(V, c(u$values, (u$values[E[, 1]] + u$values[E[, 2]]) / 2), name = u$name))
  vals <- evaluate(u, V$dof_coords)
  if (anyNA(vals)) stopf("The target mesh is not covered by the source mesh.")
  fe_function(V, vals, name = u$name)
}
