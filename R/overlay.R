# Overlay (common refinement) of two triangulations of the same domain.
#
# Every nonempty intersection of a triangle of `a` with a triangle of `b` is a
# convex polygon; it is split into a fan of triangles. On each such piece both
# piecewise polynomials are polynomials, so quadrature on the pieces
# integrates products of functions from the two meshes exactly (up to the
# quadrature degree). Clipping works in a frame translated to the lower left
# corner of both meshes, so translated meshes keep their local accuracy.
# Returns a list with the parent triangles `ea`, `eb`, the piece vertices in
# that local frame as an n by 6 matrix (x1, y1, x2, y2, x3, y3), the frame
# `origin` and the parent areas `area_a`, `area_b`.

clip_convex <- function(px, py, cx, cy) {
  # Clip polygon (px, py) by the counterclockwise triangle (cx, cy)
  # (Sutherland-Hodgman). Points on the edges are kept.
  for (i in 1:3) {
    if (!length(px)) break
    j <- i %% 3L + 1L
    ex <- cx[j] - cx[i]; ey <- cy[j] - cy[i]
    side <- ex * (py - cy[i]) - ey * (px - cx[i])
    tol <- 1e-13 * (abs(ex) + abs(ey)) * max(abs(px - cx[i]) + abs(py - cy[i]), 1e-300)
    inside <- side >= -tol
    n <- length(px)
    ox <- numeric(0); oy <- numeric(0)
    for (k in seq_len(n)) {
      l <- if (k == n) 1L else k + 1L
      if (inside[k]) { ox <- c(ox, px[k]); oy <- c(oy, py[k]) }
      if (inside[k] != inside[l]) {
        s <- side[k] / (side[k] - side[l])
        ox <- c(ox, px[k] + s * (px[l] - px[k])); oy <- c(oy, py[k] + s * (py[l] - py[k]))
      }
    }
    px <- ox; py <- oy
  }
  list(x = px, y = py)
}

mesh_overlay <- function(a, b) {
  o <- c(min(a$nodes[, 1], b$nodes[, 1]), min(a$nodes[, 2], b$nodes[, 2]))
  Pa <- shift_nodes(a$nodes, o); Ta <- a$triangles; Pb <- shift_nodes(b$nodes, o); Tb <- b$triangles
  ma <- nrow(Ta); mb <- nrow(Tb)
  bx <- matrix(Pb[Tb, 1], mb); by <- matrix(Pb[Tb, 2], mb)
  bxmin <- pmin(bx[, 1], bx[, 2], bx[, 3]); bxmax <- pmax(bx[, 1], bx[, 2], bx[, 3])
  bymin <- pmin(by[, 1], by[, 2], by[, 3]); bymax <- pmax(by[, 1], by[, 2], by[, 3])
  # Bin the triangles of b by their bounding boxes.
  xr <- range(c(Pa[, 1], Pb[, 1])); yr <- range(c(Pa[, 2], Pb[, 2]))
  nbin <- max(1L, ceiling(sqrt(mb / 2)))
  wx <- max(diff(xr), 1e-300) / nbin; wy <- max(diff(yr), 1e-300) / nbin
  bin <- function(v, lo, w) pmin(nbin - 1L, pmax(0L, floor((v - lo) / w)))
  i0 <- bin(bxmin, xr[1], wx); i1 <- bin(bxmax, xr[1], wx)
  j0 <- bin(bymin, yr[1], wy); j1 <- bin(bymax, yr[1], wy)
  cells <- vector("list", nbin * nbin)
  for (t in seq_len(mb)) {
    for (ii in i0[t]:i1[t]) for (jj in j0[t]:j1[t]) {
      k <- ii + jj * nbin + 1L
      cells[[k]] <- c(cells[[k]], t)
    }
  }
  area_b <- abs((bx[, 2] - bx[, 1]) * (by[, 3] - by[, 1]) - (bx[, 3] - bx[, 1]) * (by[, 2] - by[, 1])) / 2
  ax <- matrix(Pa[Ta, 1], ma); ay <- matrix(Pa[Ta, 2], ma)
  # One record (parents and vertices) per piece, bound once at the end;
  # growing the parent vectors piece by piece copied them every time.
  out_v <- list()
  for (e in seq_len(ma)) {
    xa <- ax[e, ]; ya <- ay[e, ]
    cand <- unique(unlist(cells[outer(bin(min(xa), xr[1], wx):bin(max(xa), xr[1], wx),
                                      bin(min(ya), yr[1], wy):bin(max(ya), yr[1], wy) * nbin, `+`) + 1L]))
    cand <- cand[bxmax[cand] >= min(xa) & bxmin[cand] <= max(xa) &
                 bymax[cand] >= min(ya) & bymin[cand] <= max(ya)]
    area_a <- abs((xa[2] - xa[1]) * (ya[3] - ya[1]) - (xa[3] - xa[1]) * (ya[2] - ya[1])) / 2
    for (t in cand) {
      pg <- clip_convex(xa, ya, bx[t, ], by[t, ])
      nv <- length(pg$x)
      if (nv < 3L) next
      # Fan triangulation; drop only slivers below rounding relative to the
      # smaller parent, so a small real cell inside a large one is kept.
      small <- 1e-12 * min(area_a, area_b[t])
      for (k in 2:(nv - 1L)) {
        v <- c(pg$x[1], pg$y[1], pg$x[k], pg$y[k], pg$x[k + 1L], pg$y[k + 1L])
        ar <- ((v[3] - v[1]) * (v[6] - v[2]) - (v[5] - v[1]) * (v[4] - v[2])) / 2
        if (abs(ar) <= small) next
        out_v[[length(out_v) + 1L]] <- c(e, t, v)
      }
    }
  }
  R <- if (length(out_v)) do.call(rbind, out_v) else matrix(numeric(0), 0, 8)
  out_ea <- as.integer(R[, 1]); out_eb <- as.integer(R[, 2]); V <- R[, 3:8, drop = FALSE]
  aa <- abs((ax[, 2] - ax[, 1]) * (ay[, 3] - ay[, 1]) - (ax[, 3] - ax[, 1]) * (ay[, 2] - ay[, 1])) / 2
  list(ea = out_ea, eb = out_eb, vertices = V, origin = o, area_a = aa, area_b = area_b)
}

shift_nodes <- function(P, o) cbind(P[, 1] - o[1], P[, 2] - o[2])

# A finite-element function whose mesh is expressed in the overlay frame.
shift_function <- function(u, o) {
  u$space$mesh$nodes <- shift_nodes(u$space$mesh$nodes, o)
  u
}

# The two meshes cover the same domain exactly when every triangle of each
# is covered by the overlay pieces; checking triangle by triangle, instead of
# total areas, also catches a small missing component or hole.
check_overlay_cover <- function(ov, what = "the error is defined only on a common domain") {
  V <- ov$vertices
  ar <- abs((V[, 3] - V[, 1]) * (V[, 6] - V[, 2]) - (V[, 5] - V[, 1]) * (V[, 4] - V[, 2])) / 2
  for (side in c("a", "b")) {
    area <- ov[[paste0("area_", side)]]
    cov <- tabulate_sum(ov[[paste0("e", side)]], ar, length(area))
    bad <- which(abs(cov - area) > 1e-8 * area)
    if (length(bad)) {
      stopf(paste0("The two meshes do not cover the same domain: triangle %d of the %s mesh is covered to ",
                   "%s of its area; %s."), bad[1], if (side == "a") "first" else "second",
            fmt_num(cov[bad[1]] / area[bad[1]]), what)
    }
  }
  ar
}

# Squared error contributions of u - ref on the overlay, summed per triangle
# of u's mesh.
overlay_error <- function(u, ref, type, diffusion, quad_degree) {
  ov <- mesh_overlay(u$space$mesh, ref$space$mesh)
  ma <- nrow(u$space$mesh$triangles)
  V <- ov$vertices
  ar <- check_overlay_cover(ov)
  o <- ov$origin
  u <- shift_function(u, o); ref <- shift_function(ref, o)
  q <- quadrature_triangle(quad_degree)
  e2 <- numeric(length(ar))
  for (k in seq_along(q$weights)) {
    l2 <- q$points[k, 1]; l3 <- q$points[k, 2]; l1 <- 1 - l2 - l3
    x <- l1 * V[, 1] + l2 * V[, 3] + l3 * V[, 5]
    y <- l1 * V[, 2] + l2 * V[, 4] + l3 * V[, 6]
    gu <- grad_in_elements(u, ov$ea, x, y)
    gr <- grad_in_elements(ref, ov$eb, x, y)
    val <- 0
    if (type %in% c("L2", "H1")) val <- val + (gu$value - gr$value)^2
    if (type %in% c("H1_semi", "H1", "energy")) {
      ex <- gu$dx - gr$dx; ey <- gu$dy - gr$dy
      d2 <- if (type != "energy") ex^2 + ey^2 else energy_density(diffusion, ex, ey, x + o[1], y + o[2])
      val <- val + d2
    }
    e2 <- e2 + 2 * q$weights[k] * ar * val
  }
  out <- numeric(ma)
  s <- rowsum(e2, ov$ea, reorder = TRUE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' Conservative L2 projection
#'
#' Computes the L2 projection `P f` of `f` onto `space`: the finite-element
#' function with `int (P f) v = int f v` for every `v` in the space. Because
#' constants belong to the space, the projection preserves the integral of
#' `f` over the whole domain, which interpolation does not. It does not
#' preserve integrals over individual elements or unions of them: local
#' conservation would need a discontinuous space or explicit constraints.
#'
#' `f` may be a number, a `function(x, y)`, or an [fe_function()] on any mesh
#' of the same domain. For a finite-element function on another mesh the
#' load integrals are computed exactly on the overlay of the two meshes
#' (pointwise on the target mesh when it refines the source mesh), so the
#' transfer is conservative up to rounding.
#'
#' @param f Function to project.
#' @param space Target [fe_space()].
#' @param quad_degree Quadrature exactness degree for the load integrals.
#' @return An [fe_function()] on `space`.
#' @examples
#' a <- interpolate(function(x, y) sin(5 * x) * y, fe_space(pde_mesh_rectangle(nx = 12)))
#' V <- fe_space(pde_mesh_rectangle(nx = 5, pattern = "left"))
#' p <- l2_project(a, V)
#' c(sum(assemble_matrix(a$space, "mass") %*% a$values),
#'   sum(assemble_matrix(V, "mass") %*% p$values))
#' @export
l2_project <- function(f, space, quad_degree = NULL) {
  stopifnot(inherits(space, "fe_space"))
  # The mass matrix is always integrated exactly (degree 2p products);
  # `quad_degree` applies to the load only.
  M <- assemble_matrix(space, "mass", quad_degree = max(2L * space$degree, if (is.null(quad_degree)) 0L else quad_degree))
  b <- if (inherits(f, "fe_function") && !same_space(f$space, space) &&
           is.null(nested_parent(space$mesh, f$space$mesh))) {
    if (is.null(quad_degree)) quad_degree <- f$space$degree + space$degree + 2L
    overlay_load(f, space, quad_degree)
  } else {
    assemble_vector(space, f, quad_degree = quad_degree)
  }
  s <- solve_linear_system(M, b)
  warn_solver(s$info, "projection solve")
  fe_function(space, s$x, name = if (inherits(f, "fe_function")) f$name else NULL)
}

# int f phi_i over the overlay of f's mesh and the target mesh.
overlay_load <- function(f, space, quad_degree) {
  ov <- mesh_overlay(space$mesh, f$space$mesh)
  V <- ov$vertices
  ar <- check_overlay_cover(ov, "the projection is defined only between meshes of one domain")
  f <- shift_function(f, ov$origin)
  P <- shift_nodes(space$mesh$nodes, ov$origin); T <- space$mesh$triangles[ov$ea, , drop = FALSE]
  x1 <- P[T[, 1], 1]; y1 <- P[T[, 1], 2]
  b11 <- P[T[, 2], 1] - x1; b12 <- P[T[, 3], 1] - x1
  b21 <- P[T[, 2], 2] - y1; b22 <- P[T[, 3], 2] - y1
  det <- b11 * b22 - b12 * b21
  q <- quadrature_triangle(quad_degree)
  nl <- space$element$ndof
  B <- matrix(0, length(ar), nl)
  for (k in seq_along(q$weights)) {
    l2 <- q$points[k, 1]; l3 <- q$points[k, 2]; l1 <- 1 - l2 - l3
    x <- l1 * V[, 1] + l2 * V[, 3] + l3 * V[, 5]
    y <- l1 * V[, 2] + l2 * V[, 4] + l3 * V[, 6]
    fv <- grad_in_elements(f, ov$eb, x, y)$value
    dx <- x - x1; dy <- y - y1
    xi <- (b22 * dx - b12 * dy) / det; eta <- (-b21 * dx + b11 * dy) / det
    B <- B + (2 * q$weights[k] * ar * fv) * space$element$basis(xi, eta)
  }
  out <- numeric(space$ndof)
  s <- rowsum(as.vector(B), as.vector(space$cell_dofs[ov$ea, , drop = FALSE]), reorder = TRUE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}
