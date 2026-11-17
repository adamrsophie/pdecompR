# Element geometry -----------------------------------------------------------

# Precompute, for every element and quadrature point, the physical
# coordinates, basis values and physical basis gradients.
element_geometry <- function(space, quad_degree = NULL) {
  el <- space$element
  if (is.null(quad_degree)) quad_degree <- 2L * el$degree + 2L
  q <- quadrature_triangle(quad_degree)
  g <- affine_geometry(space$mesh, q)
  xi <- q$points[, 1]; eta <- q$points[, 2]
  phi <- el$basis(xi, eta)
  G <- el$grad(xi, eta)
  m <- length(g$det)
  nq <- length(q$weights)
  gx <- gy <- vector("list", nq)
  # P1 gradients are constant on each triangle: compute them once and share
  # the same matrix across quadrature points.
  same <- all(G$xi == rep(G$xi[1, ], each = nq)) && all(G$eta == rep(G$eta[1, ], each = nq))
  for (k in seq_len(nq)) {
    if (same && k > 1L) { gx[[k]] <- gx[[1]]; gy[[k]] <- gy[[1]]; next }
    gxi <- matrix(G$xi[k, ], m, el$ndof, byrow = TRUE)
    geta <- matrix(G$eta[k, ], m, el$ndof, byrow = TRUE)
    gx[[k]] <- (g$b22 * gxi - g$b21 * geta) / g$det
    gy[[k]] <- (-g$b12 * gxi + g$b11 * geta) / g$det
  }
  list(space = space, quad = q, w = q$weights, det = g$det, x = g$x, y = g$y,
       phi = phi, gx = gx, gy = gy, nq = nq, nloc = el$ndof, m = m,
       b11 = g$b11, b12 = g$b12, b21 = g$b21, b22 = g$b22, gxi = G$xi, geta = G$eta)
}

# Evaluate a scalar coefficient at all quadrature points (m by nq matrix).
coef_qp <- function(coef, geo, t = NULL, name = "coefficient") {
  if (is.null(coef)) return(NULL)
  if (inherits(coef, "pde_region_coefficient")) {
    m0 <- attr(coef, "mesh")
    if (identical(m0$nodes, geo$space$mesh$nodes) && identical(m0$triangles, geo$space$mesh$triangles)) {
      return(matrix(attr(coef, "element_values"), geo$m, geo$nq))
    }
  }
  if (inherits(coef, "fe_function")) {
    if (same_space(coef$space, geo$space)) return(fe_qp(coef, geo)$value)
    # Another degree or another mesh: evaluate the function at the points.
    v <- evaluate(coef, as.vector(geo$x), as.vector(geo$y))
    if (anyNA(v)) stopf("`%s` (a finite-element function) does not cover the whole domain.", name)
    return(matrix(v, geo$m, geo$nq))
  }
  matrix(eval_spatial(coef, as.vector(geo$x), as.vector(geo$y), t = t, name = name),
         geo$m, geo$nq)
}

# Evaluate a vector coefficient at all quadrature points; returns list(x, y).
vec_qp <- function(b, geo, t = NULL, name = "advection") {
  if (is.null(b)) return(NULL)
  if (is.numeric(b)) {
    if (length(b) != 2L) stopf("`%s` must be a vector of length 2 or a function.", name)
    return(list(x = matrix(b[1], geo$m, geo$nq), y = matrix(b[2], geo$m, geo$nq)))
  }
  if (!is.function(b)) stopf("`%s` must be a vector of length 2 or a function.", name)
  xx <- as.vector(geo$x); yy <- as.vector(geo$y)
  val <- if (!is.null(t) && is_time_function(b)) b(xx, yy, t) else b(xx, yy)
  val <- as.matrix(val)
  if (nrow(val) != length(xx) || ncol(val) != 2L) {
    stopf("`%s` must return a matrix with two columns and one row per point.", name)
  }
  list(x = matrix(val[, 1], geo$m, geo$nq), y = matrix(val[, 2], geo$m, geo$nq))
}

same_space <- function(a, b) {
  a$degree == b$degree && identical(a$mesh$nodes, b$mesh$nodes) &&
    identical(a$mesh$triangles, b$mesh$triangles)
}

# Values and gradients of a finite-element function at quadrature points.
fe_qp <- function(u, geo) {
  vals <- if (inherits(u, "fe_function")) u$values else u
  U <- matrix(vals[geo$space$cell_dofs], geo$m)
  value <- U %*% t(geo$phi)
  dx <- dy <- matrix(0, geo$m, geo$nq)
  for (k in seq_len(geo$nq)) {
    dx[, k] <- rowSums(geo$gx[[k]] * U)
    dy[, k] <- rowSums(geo$gy[[k]] * U)
  }
  list(value = value, dx = dx, dy = dy)
}

# Local bilinear kernel --------------------------------------------------------
#
# For every element computes
#   A[i, j] = sum_q w_q det [ mass phi_i phi_j
#                             + phi_i (adv . grad phi_j)
#                             + (advt . grad phi_i) phi_j
#                             + grad phi_i^T D grad phi_j ]
# with i the test index (row) and j the trial index (column). `diff` is an
# m by nq scalar matrix or a list(xx, xy, yx, yy) of such matrices.
local_bilinear <- function(geo, mass = NULL, adv = NULL, advt = NULL, diff = NULL) {
  m <- geo$m; nl <- geo$nloc
  A <- array(0, c(m, nl, nl))
  # Terms whose coefficients are constant on every element are assembled
  # from reference-element matrices: on affine triangles the quadrature sum
  # factorizes exactly into (per-element metric) x (reference integral), so
  # the result equals the quadrature result up to rounding at a fraction of
  # the cost. Other terms go through the quadrature loop below.
  pec <- function(a) is.matrix(a) && isTRUE(all(a == a[, 1]))
  if (!is.null(geo$gxi)) {
    w <- geo$w; det <- geo$det; ph <- geo$phi
    ref <- function(a, b) crossprod(a, w * b)
    G11 <- geo$b22 / det; G12 <- -geo$b21 / det; G21 <- -geo$b12 / det; G22 <- geo$b11 / det
    add <- function(cf, S) A <<- A + array(cf %o% S, c(m, nl, nl))
    if (!is.null(mass) && pec(mass)) { add(det * mass[, 1], ref(ph, ph)); mass <- NULL }
    if (!is.null(diff) && (if (is.list(diff)) all(vapply(diff, pec, logical(1))) else pec(diff))) {
      d <- if (is.list(diff)) lapply(diff, function(a) a[, 1]) else list(xx = diff[, 1], xy = 0, yx = 0, yy = diff[, 1])
      DG11 <- d$xx * G11 + d$xy * G21; DG12 <- d$xx * G12 + d$xy * G22
      DG21 <- d$yx * G11 + d$yy * G21; DG22 <- d$yx * G12 + d$yy * G22
      gxi <- geo$gxi; geta <- geo$geta
      add(det * (G11 * DG11 + G21 * DG21), ref(gxi, gxi))
      add(det * (G11 * DG12 + G21 * DG22), ref(gxi, geta))
      add(det * (G12 * DG11 + G22 * DG21), ref(geta, gxi))
      add(det * (G12 * DG12 + G22 * DG22), ref(geta, geta))
      diff <- NULL
    }
    if (!is.null(adv) && pec(adv$x) && pec(adv$y)) {
      bx <- adv$x[, 1]; by <- adv$y[, 1]
      add(det * (G11 * bx + G21 * by), ref(ph, geo$gxi))
      add(det * (G12 * bx + G22 * by), ref(ph, geo$geta))
      adv <- NULL
    }
    if (!is.null(advt) && pec(advt$x) && pec(advt$y)) {
      bx <- advt$x[, 1]; by <- advt$y[, 1]
      add(det * (G11 * bx + G21 * by), t(ref(ph, geo$gxi)))
      add(det * (G12 * bx + G22 * by), t(ref(ph, geo$geta)))
      advt <- NULL
    }
    if (is.null(mass) && is.null(diff) && is.null(adv) && is.null(advt)) return(A)
  }
  tensor <- is.list(diff)
  for (k in seq_len(geo$nq)) {
    wd <- geo$w[k] * geo$det
    gx <- geo$gx[[k]]; gy <- geo$gy[[k]]; ph <- geo$phi[k, ]
    for (j in seq_len(nl)) {
      for (i in seq_len(nl)) {
        v <- 0
        if (!is.null(mass)) v <- v + mass[, k] * ph[i] * ph[j]
        if (!is.null(adv)) v <- v + ph[i] * (adv$x[, k] * gx[, j] + adv$y[, k] * gy[, j])
        if (!is.null(advt)) v <- v + (advt$x[, k] * gx[, i] + advt$y[, k] * gy[, i]) * ph[j]
        if (!is.null(diff)) {
          if (tensor) {
            v <- v + gx[, i] * (diff$xx[, k] * gx[, j] + diff$xy[, k] * gy[, j]) +
              gy[, i] * (diff$yx[, k] * gx[, j] + diff$yy[, k] * gy[, j])
          } else {
            v <- v + diff[, k] * (gx[, i] * gx[, j] + gy[, i] * gy[, j])
          }
        }
        A[, i, j] <- A[, i, j] + wd * v
      }
    }
  }
  A
}

# Local linear kernel: b[i] = sum_q w det (f phi_i + g . grad phi_i).
local_linear <- function(geo, f = NULL, g = NULL) {
  m <- geo$m; nl <- geo$nloc
  b <- matrix(0, m, nl)
  if (!is.null(f) && is.null(g) && isTRUE(all(f == f[, 1]))) {
    return((geo$det * f[, 1]) %o% colSums(geo$w * geo$phi))
  }
  for (k in seq_len(geo$nq)) {
    wd <- geo$w[k] * geo$det
    for (i in seq_len(nl)) {
      v <- 0
      if (!is.null(f)) v <- v + f[, k] * geo$phi[k, i]
      if (!is.null(g)) v <- v + g$x[, k] * geo$gx[[k]][, i] + g$y[, k] * geo$gy[[k]][, i]
      b[, i] <- b[, i] + wd * v
    }
  }
  b
}

global_matrix <- function(space, A) {
  D <- space$cell_dofs
  nl <- ncol(D)
  I <- D[, rep(seq_len(nl), times = nl)]
  J <- D[, rep(seq_len(nl), each = nl)]
  Matrix::sparseMatrix(i = as.vector(I), j = as.vector(J), x = as.vector(A),
                       dims = c(space$ndof, space$ndof))
}

# Assembly plan for repeated assembly on one space: the sparsity pattern
# and the map from element triplets to its nonzeros are computed once; each
# call then only sums the element values into the nonzeros (one sparse
# matrix-vector product), which matters when coefficients change every
# time step or Newton iteration.
assembly_plan <- function(space) {
  D <- space$cell_dofs; nl <- ncol(D); n <- space$ndof
  i <- as.vector(D[, rep(seq_len(nl), times = nl)])
  j <- as.vector(D[, rep(seq_len(nl), each = nl)])
  key <- (j - 1) * n + i
  u <- sort(unique(key))
  S <- Matrix::sparseMatrix(i = match(key, u), j = seq_along(key), x = 1, dims = c(length(u), length(key)))
  tmpl <- Matrix::sparseMatrix(i = (u - 1) %% n + 1, j = (u - 1) %/% n + 1, x = rep(1, length(u)), dims = c(n, n))
  function(A) {
    out <- tmpl
    out@x <- as.numeric(S %*% as.vector(A))
    out
  }
}

global_vector <- function(space, b) {
  out <- numeric(space$ndof)
  s <- rowsum(as.vector(b), as.vector(space$cell_dofs), reorder = TRUE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' Finite-element assembly
#'
#' `assemble_matrix()` assembles a global sparse matrix from element
#' contributions using sparse triplets, which are summed into a compressed
#' column `dgCMatrix`. Element integrals use triangle quadrature exact to
#' `quad_degree` (default `2 p + 2` for degree `p`), which is exact for
#' constant and affine coefficients.
#'
#' Matrix types (`i` is the test function, `j` the trial function):
#'
#' * `"stiffness"`: `K_ij = int k grad(phi_j) . grad(phi_i)`; `coef` is the
#'   diffusion `k` (number, function or 2 by 2 matrix for a constant tensor).
#' * `"mass"`: `M_ij = int c phi_j phi_i`; `coef` defaults to 1.
#' * `"advection"`: `C_ij = int (b . grad(phi_j)) phi_i`; `coef` is the
#'   velocity `b` (length-2 vector or `function(x, y)` returning a two
#'   column matrix).
#'
#' `assemble_vector()` assembles `F_i = int f phi_i`.
#' `assemble_boundary_matrix()` and `assemble_boundary_vector()` assemble
#' `int_Gamma alpha phi_j phi_i ds` and `int_Gamma g phi_i ds` over boundary
#' edges with the given markers, using Gauss-Legendre rules on each edge.
#'
#' @param space An [fe_space()].
#' @param type Matrix type.
#' @param coef Coefficient, see Details.
#' @param f,g,alpha Scalar coefficients (numbers or `function(x, y)`, or
#'   `function(x, y, t)` together with `t`).
#' @param marker Boundary markers.
#' @param t Optional time passed to time-dependent coefficients.
#' @param quad_degree Quadrature exactness degree.
#' @return A sparse `dgCMatrix` or a numeric vector.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 4))
#' K <- assemble_matrix(V, "stiffness")
#' M <- assemble_matrix(V, "mass")
#' sum(M)                                  # area of the domain
#' max(abs(K %*% rep(1, V$ndof)))          # constants are in the kernel
#' @export
assemble_matrix <- function(space, type = c("stiffness", "mass", "advection"),
                            coef = NULL, t = NULL, quad_degree = NULL) {
  type <- match.arg(type)
  geo <- element_geometry(space, quad_degree)
  A <- switch(type,
    stiffness = local_bilinear(geo, diff = diffusion_qp(if (is.null(coef)) 1 else coef, geo, t)),
    mass = local_bilinear(geo, mass = coef_qp(if (is.null(coef)) 1 else coef, geo, t, "coef")),
    advection = {
      if (is.null(coef)) stopf("`coef` (the velocity) is required for advection.")
      local_bilinear(geo, adv = vec_qp(coef, geo, t))
    })
  global_matrix(space, A)
}

#' @rdname assemble_matrix
#' @export
assemble_vector <- function(space, f, t = NULL, quad_degree = NULL) {
  geo <- element_geometry(space, quad_degree)
  global_vector(space, local_linear(geo, f = coef_qp(f, geo, t, "f")))
}

#' @rdname assemble_matrix
#' @export
assemble_boundary_matrix <- function(space, alpha = 1, marker = "boundary",
                                     t = NULL, quad_degree = NULL) {
  bg <- boundary_geometry(space, marker, quad_degree)
  if (is.null(bg)) return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                                dims = c(space$ndof, space$ndof)))
  a <- matrix(eval_spatial(alpha, as.vector(bg$x), as.vector(bg$y), t, "alpha"), bg$ne)
  nl <- bg$nloc
  A <- array(0, c(bg$ne, nl, nl))
  for (k in seq_along(bg$w)) {
    wl <- bg$w[k] * bg$len * a[, k]
    for (j in seq_len(nl)) for (i in seq_len(nl)) {
      A[, i, j] <- A[, i, j] + wl * bg$phi[[k]][, i] * bg$phi[[k]][, j]
    }
  }
  D <- bg$dofs
  Matrix::sparseMatrix(i = as.vector(D[, rep(seq_len(nl), times = nl)]),
                       j = as.vector(D[, rep(seq_len(nl), each = nl)]),
                       x = as.vector(A), dims = c(space$ndof, space$ndof))
}

#' @rdname assemble_matrix
#' @export
assemble_boundary_vector <- function(space, g, marker = "boundary", t = NULL,
                                     quad_degree = NULL) {
  bg <- boundary_geometry(space, marker, quad_degree)
  out <- numeric(space$ndof)
  if (is.null(bg)) return(out)
  gv <- matrix(eval_spatial(g, as.vector(bg$x), as.vector(bg$y), t, "g"), bg$ne)
  b <- matrix(0, bg$ne, bg$nloc)
  for (k in seq_along(bg$w)) {
    b <- b + (bg$w[k] * bg$len * gv[, k]) * bg$phi[[k]]
  }
  s <- rowsum(as.vector(b), as.vector(bg$dofs), reorder = TRUE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' Spatially varying diffusion tensors
#'
#' Builds a symmetric diffusion tensor `K = [[xx, xy], [xy, yy]]` whose
#' entries are numbers or functions of `(x, y)` (or `(x, y, t)`), for
#' anisotropic media whose principal directions or strengths vary in space.
#' Use it wherever a diffusion coefficient is accepted. Positive
#' definiteness is checked at every quadrature point during assembly.
#'
#' @param xx,xy,yy Tensor entries.
#' @return An object of class `pde_tensor`.
#' @examples
#' # Diffusion aligned with circles around the origin
#' K <- tensor_coefficient(function(x, y) 1 + y^2, function(x, y) -x * y,
#'                         function(x, y) 1 + x^2)
#' p <- elliptic_problem(pde_mesh_rectangle(nx = 8), diffusion = K, source = 1,
#'                       bcs = dirichlet_bc(0))
#' max(solve_pde(p)$u$values)
#' @export
tensor_coefficient <- function(xx, xy = 0, yy = xx) {
  for (nm in c("xx", "xy", "yy")) check_coefficient(get(nm), paste0("tensor entry `", nm, "`"))
  structure(list(xx = xx, xy = xy, yy = yy), class = "pde_tensor")
}

#' @export
print.pde_tensor <- function(x, ...) {
  d <- function(v) if (is.function(v)) "function" else fmt_num(v)
  cat("<pde_tensor> xx = ", d(x$xx), ", xy = ", d(x$xy), ", yy = ", d(x$yy), "\n", sep = "")
  invisible(x)
}

# Entries of a constant 2 by 2 matrix or a pde_tensor at points, or NULL for
# a scalar coefficient.
tensor_at <- function(k, x, y, t = NULL) {
  if (inherits(k, "pde_tensor")) {
    e <- lapply(k[c("xx", "xy", "yy")], function(v) eval_spatial(v, x, y, t = t, name = "diffusion tensor"))
    return(list(xx = e$xx, xy = e$xy, yx = e$xy, yy = e$yy))
  }
  if (is.matrix(k) && is.numeric(k)) {
    n <- length(x)
    return(list(xx = rep(k[1, 1], n), xy = rep(k[1, 2], n), yx = rep(k[2, 1], n), yy = rep(k[2, 2], n)))
  }
  NULL
}

# K-weighted squared norm of (ex, ey) for scalar or tensor diffusion.
energy_density <- function(k, ex, ey, x, y, kq = NULL) {
  K <- tensor_at(k, x, y)
  if (is.null(K)) return((ex^2 + ey^2) * (if (is.null(kq)) eval_spatial(k, x, y, name = "diffusion") else kq))
  K$xx * ex^2 + (K$xy + K$yx) * ex * ey + K$yy * ey^2
}

diffusion_qp <- function(k, geo, t = NULL) {
  if (inherits(k, "pde_tensor")) {
    K <- tensor_at(k, as.vector(geo$x), as.vector(geo$y), t)
    if (any(!is.finite(unlist(K)))) stopf("The diffusion tensor must be finite.")
    if (any(K$xx <= 0 | K$xx * K$yy - K$xy^2 <= 0)) {
      stopf("The diffusion tensor must be positive definite at every quadrature point.")
    }
    return(lapply(K, matrix, geo$m, geo$nq))
  }
  if (is.matrix(k) && is.numeric(k)) {
    if (!all(dim(k) == 2L)) stopf("A tensor diffusion coefficient must be a 2 by 2 matrix.")
    if (any(!is.finite(k))) stopf("The diffusion tensor must be finite.")
    if (abs(k[1, 2] - k[2, 1]) > 1e-12 * max(abs(k))) stopf("The diffusion tensor must be symmetric.")
    ev <- eigen((k + t(k)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) stopf("The diffusion tensor must be positive definite.")
    z <- matrix(1, geo$m, geo$nq)
    return(list(xx = k[1, 1] * z, xy = k[1, 2] * z, yx = k[2, 1] * z, yy = k[2, 2] * z))
  }
  coef_qp(k, geo, t, "diffusion")
}

# Boundary edge geometry: quadrature points, lengths, normals and the values
# of the parent element's basis functions at the edge quadrature points.
boundary_geometry <- function(space, marker, quad_degree = NULL) {
  mesh <- space$mesh
  rows <- boundary_rows(mesh, marker)
  if (!length(rows)) return(NULL)
  bd <- mesh$boundary[rows, , drop = FALSE]
  el <- space$element
  if (is.null(quad_degree)) quad_degree <- 2L * el$degree + 2L
  gq <- quadrature_interval(ceiling((quad_degree + 1) / 2))
  s <- gq$points[, 1]
  P <- mesh$nodes
  x1 <- P[bd$n1, 1]; y1 <- P[bd$n1, 2]; x2 <- P[bd$n2, 1]; y2 <- P[bd$n2, 2]
  len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  start <- rbind(c(0, 0), c(1, 0), c(0, 1)); end <- rbind(c(1, 0), c(0, 1), c(0, 0))
  phi <- vector("list", length(s))
  for (k in seq_along(s)) {
    xi <- (1 - s[k]) * start[bd$local, 1] + s[k] * end[bd$local, 1]
    eta <- (1 - s[k]) * start[bd$local, 2] + s[k] * end[bd$local, 2]
    phi[[k]] <- el$basis(xi, eta)
  }
  list(ne = nrow(bd), nloc = el$ndof, w = gq$weights, len = len,
       x = x1 + outer(x2 - x1, s), y = y1 + outer(y2 - y1, s),
       nx = (y2 - y1) / len, ny = -(x2 - x1) / len,
       phi = phi, dofs = space$cell_dofs[bd$element, , drop = FALSE],
       element = bd$element, s = s, rows = rows)
}
