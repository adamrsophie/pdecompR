#' Triangular meshes
#'
#' `pde_mesh()` builds and validates a two-dimensional triangular mesh from
#' node coordinates and element connectivity.
#'
#' @section Conventions:
#' * Indices are one-based, as everywhere in R. Node `i` is row `i` of
#'   `nodes`; element `e` is row `e` of `triangles`.
#' * Triangles are stored counterclockwise (positive signed area).
#' * The local edges of a triangle `(v1, v2, v3)` are edge 1 `(v1, v2)`,
#'   edge 2 `(v2, v3)` and edge 3 `(v3, v1)`.
#' * Boundary edges are oriented so that the domain lies on their left;
#'   their outward normal is the edge direction rotated clockwise.
#' * Each boundary edge carries exactly one marker (a character label). The
#'   label `"boundary"` always selects all boundary edges.
#'
#' @section Validation:
#' The constructor rejects non-finite coordinates, invalid or repeated node
#' indices, duplicate elements, degenerate (zero-area) triangles, unused nodes,
#' edges shared by more than two triangles, and boundary definitions that do
#' not match the topological boundary. Clockwise triangles are rejected unless
#' `fix_orientation = TRUE`, in which case they are reordered and the number of
#' reordered elements is reported.
#'
#' @param nodes Numeric matrix with two columns (x, y).
#' @param triangles Integer matrix with three columns of node indices.
#' @param boundary_edges Optional two-column matrix of node pairs on the
#'   boundary. When given, it must list every boundary edge exactly once.
#' @param boundary_markers Character vector with one marker per row of
#'   `boundary_edges`. Ignored when `boundary_edges` is `NULL`; all boundary
#'   edges are then marked `"boundary"`.
#' @param fix_orientation Reorder clockwise triangles instead of failing.
#' @param regions Optional character vector with one region (subdomain) name
#'   per triangle, for piecewise coefficients ([region_coefficient()]).
#'   Default: every triangle belongs to region `"domain"`.
#' @return An object of class `pde_mesh` with components `nodes`, `triangles`,
#'   `edges`, `element_edges`, `edge_elements`, `boundary` (a data frame of
#'   boundary edges) and `info`.
#' @examples
#' nodes <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' tris <- rbind(c(1, 2, 3), c(1, 3, 4))
#' mesh <- pde_mesh(nodes, tris)
#' mesh
#' @export
pde_mesh <- function(nodes, triangles, boundary_edges = NULL,
                     boundary_markers = NULL, fix_orientation = FALSE, regions = NULL) {
  nodes <- as.matrix(nodes)
  if (!is.numeric(nodes) || ncol(nodes) != 2L) {
    stopf("`nodes` must be a numeric matrix with two columns.")
  }
  storage.mode(nodes) <- "double"
  dimnames(nodes) <- NULL
  if (any(!is.finite(nodes))) stopf("`nodes` contains non-finite coordinates.")
  tris <- as.matrix(triangles)
  if (!is.numeric(tris) || ncol(tris) != 3L || nrow(tris) < 1L) {
    stopf("`triangles` must be a numeric matrix with three columns.")
  }
  if (any(!is.finite(tris)) || any(tris != round(tris))) {
    stopf("`triangles` must contain whole-number node indices.")
  }
  storage.mode(tris) <- "integer"
  dimnames(tris) <- NULL
  n <- nrow(nodes)
  if (any(tris < 1L | tris > n)) {
    bad <- which(rowSums(tris < 1L | tris > n) > 0)[1]
    stopf("Element %d refers to a node outside 1..%d.", bad, n)
  }
  rep_node <- tris[, 1] == tris[, 2] | tris[, 2] == tris[, 3] | tris[, 1] == tris[, 3]
  if (any(rep_node)) {
    stopf("Element %d repeats a node index.", which(rep_node)[1])
  }
  lo <- pmin(tris[, 1], tris[, 2], tris[, 3]); hi <- pmax(tris[, 1], tris[, 2], tris[, 3])
  dup <- duplicated(cbind(lo, tris[, 1] + tris[, 2] + tris[, 3] - lo - hi, hi))
  if (any(dup)) stopf("Element %d duplicates an earlier element.", which(dup)[1])
  used <- tabulate(tris, nbins = n)
  if (any(used == 0L)) {
    stopf("%d node(s) are not used by any element (first: %d).",
          sum(used == 0L), which(used == 0L)[1])
  }
  det <- tri_det(nodes, tris)
  e2 <- function(a, b) rowSums((nodes[tris[, a], , drop = FALSE] -
                                  nodes[tris[, b], , drop = FALSE])^2)
  lmax2 <- pmax(e2(1, 2), e2(2, 3), e2(3, 1))
  degenerate <- abs(det) <= 1e-12 * lmax2
  if (any(degenerate)) {
    stopf("Element %d has zero area (degenerate triangle).", which(degenerate)[1])
  }
  neg <- det < 0
  n_fixed <- 0L
  if (any(neg)) {
    if (!fix_orientation) {
      stopf(paste0("%d element(s) are clockwise (first: %d). ",
                   "Use `fix_orientation = TRUE` to reorder them."),
            sum(neg), which(neg)[1])
    }
    tris[neg, c(2, 3)] <- tris[neg, c(3, 2)]
    n_fixed <- sum(neg)
  }
  mesh <- build_topology(nodes, tris)
  mesh$boundary <- assign_boundary(mesh, boundary_edges, boundary_markers)
  if (is.null(regions)) regions <- "domain"
  if (length(regions) == 1L) regions <- rep(regions, nrow(tris))
  if (!is.character(regions) || length(regions) != nrow(tris) || anyNA(regions) || any(!nzchar(regions))) {
    stopf("`regions` must be one nonempty name or one name per triangle.")
  }
  mesh$regions <- regions
  mesh$info <- list(reoriented = n_fixed, generator = "user")
  mesh
}

#' Subdomains (regions) of a mesh
#'
#' Every triangle of a mesh belongs to a named region (subdomain); by default
#' all triangles belong to `"domain"`. `mesh_regions()` returns the region of
#' each triangle, `set_region()` assigns a region to the triangles whose
#' centroids satisfy a predicate, and [read_gmsh()] takes regions from
#' physical surfaces. Refinement keeps the region of the parent triangle.
#' Regions carry piecewise coefficients through [region_coefficient()].
#'
#' @param mesh A [pde_mesh()].
#' @param region Region name.
#' @param which A function `function(x, y)` returning `TRUE` for the
#'   centroids of the triangles to assign.
#' @return `mesh_regions()`: a character vector with one entry per triangle.
#'   `set_region()`: the updated mesh.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 4)
#' mesh <- set_region(mesh, "inclusion", function(x, y) (x - 0.5)^2 + (y - 0.5)^2 < 0.1)
#' table(mesh_regions(mesh))
#' @export
mesh_regions <- function(mesh) {
  stopifnot(inherits(mesh, "pde_mesh"))
  if (is.null(mesh$regions)) rep("domain", nrow(mesh$triangles)) else mesh$regions
}

#' @rdname mesh_regions
#' @export
set_region <- function(mesh, region, which) {
  stopifnot(inherits(mesh, "pde_mesh"))
  check_name(region, "region")
  if (!is.function(which)) stopf("`which` must be a function(x, y).")
  P <- mesh$nodes; T <- mesh$triangles
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
  sel <- check_selection(which(cx, cy), nrow(T), "triangle")
  r <- mesh_regions(mesh)
  r[sel] <- region
  mesh$regions <- r
  mesh
}

# Give the triangles of a nested (refined) mesh the regions of their parents.
inherit_regions <- function(child, parent) {
  r <- mesh_regions(parent)
  if (length(unique(r)) <= 1L) {
    child$regions <- rep(r[1], nrow(child$triangles))
    return(child)
  }
  P <- child$nodes; T <- child$triangles
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
  child$regions <- r[locate_points(parent, cx, cy)$element]
  child
}

# Edge structures for a valid, counterclockwise triangulation.
build_topology <- function(nodes, tris) {
  n <- nrow(nodes); m <- nrow(tris)
  a <- c(tris[, 1], tris[, 2], tris[, 3])
  b <- c(tris[, 2], tris[, 3], tris[, 1])
  key <- edge_key(a, b, n)
  ukey <- unique(key)
  id <- match(key, ukey)
  ne <- length(ukey)
  first <- match(seq_len(ne), id)
  edges <- cbind(a[first], b[first])
  element_edges <- matrix(id, nrow = m, ncol = 3)
  elem <- rep(seq_len(m), 3)
  cnt <- tabulate(id, nbins = ne)
  if (any(cnt > 2L)) {
    stopf("Edge (%d, %d) is shared by more than two elements.",
          edges[which(cnt > 2L)[1], 1], edges[which(cnt > 2L)[1], 2])
  }
  # A shared edge must be traversed in opposite directions by its two
  # elements; equal directions mean overlapping, inconsistent elements.
  ord <- order(id)
  id_s <- id[ord]
  second <- ord[c(FALSE, id_s[-1] == id_s[-length(id_s)])]
  first_occ <- first[id[second]]
  if (length(second) && any(a[second] == a[first_occ])) {
    bad <- second[a[second] == a[first_occ]][1]
    stopf("Elements %d and %d traverse a shared edge in the same direction (overlapping elements).",
          elem[first[id[bad]]], elem[bad])
  }
  edge_elements <- matrix(NA_integer_, ne, 2)
  edge_elements[, 1] <- elem[first]
  edge_elements[id[second], 2] <- elem[second]
  local <- rep(1:3, each = m)
  bnd <- which(cnt == 1L)
  bpos <- first[bnd]
  structure(list(
    nodes = nodes, triangles = tris, edges = edges,
    element_edges = element_edges, edge_elements = edge_elements,
    boundary = data.frame(n1 = a[bpos], n2 = b[bpos], edge = bnd,
                          element = elem[bpos], local = local[bpos],
                          marker = "boundary", stringsAsFactors = FALSE),
    dim = 2L
  ), class = "pde_mesh")
}

assign_boundary <- function(mesh, boundary_edges, boundary_markers) {
  bd <- mesh$boundary
  if (is.null(boundary_edges)) return(bd)
  be <- as.matrix(boundary_edges)
  if (ncol(be) != 2L) stopf("`boundary_edges` must have two columns.")
  if (is.null(boundary_markers)) boundary_markers <- rep("boundary", nrow(be))
  if (length(boundary_markers) == 1L) boundary_markers <- rep(boundary_markers, nrow(be))
  if (length(boundary_markers) != nrow(be)) {
    stopf("`boundary_markers` must have one entry per boundary edge.")
  }
  n <- nrow(mesh$nodes)
  if (!is.numeric(be) || any(!is.finite(be)) || any(be != round(be)) || any(be < 1 | be > n)) {
    stopf("`boundary_edges` must contain node indices between 1 and %d.", n)
  }
  if (!is.character(boundary_markers) || anyNA(boundary_markers) || any(!nzchar(boundary_markers))) {
    stopf("`boundary_markers` must be nonempty names.")
  }
  k_given <- edge_key(be[, 1], be[, 2], n)
  if (anyDuplicated(k_given)) stopf("`boundary_edges` lists an edge twice.")
  k_topo <- edge_key(bd$n1, bd$n2, n)
  pos <- match(k_given, k_topo)
  if (anyNA(pos)) {
    i <- which(is.na(pos))[1]
    stopf("Boundary edge (%d, %d) is not on the mesh boundary.", be[i, 1], be[i, 2])
  }
  if (length(pos) != length(k_topo)) {
    stopf("%d topological boundary edge(s) have no marker (missing boundary definition).",
          length(k_topo) - length(pos))
  }
  bd$marker[pos] <- as.character(boundary_markers)
  bd
}

#' Structured triangular mesh of a rectangle
#'
#' Splits the rectangle into `nx` by `ny` cells and each cell into two
#' triangles.
#'
#' @param x,y Endpoints of the rectangle in each direction.
#' @param nx,ny Number of cells in each direction.
#' @param pattern Diagonal pattern: `"right"` (diagonals from lower left to
#'   upper right), `"left"` (lower right to upper left) or `"alternating"`
#'   (the diagonal direction alternates from cell to cell).
#' @return A `pde_mesh` with boundary markers `"bottom"`, `"right"`, `"top"`
#'   and `"left"`. Corner nodes belong to both adjacent sides.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 4, ny = 4)
#' summary(mesh)
#' @export
pde_mesh_rectangle <- function(x = c(0, 1), y = c(0, 1), nx = 8, ny = nx,
                               pattern = c("right", "left", "alternating")) {
  pattern <- match.arg(pattern)
  check_endpoints(x, "x"); check_endpoints(y, "y")
  nx <- check_count(nx, "nx"); ny <- check_count(ny, "ny")
  xs <- x[1] + (x[2] - x[1]) * (0:nx) / nx
  ys <- y[1] + (y[2] - y[1]) * (0:ny) / ny
  nodes <- cbind(rep(xs, ny + 1), rep(ys, each = nx + 1))
  id <- function(i, j) i + j * (nx + 1L) + 1L
  ci <- rep(0:(nx - 1L), ny); cj <- rep(0:(ny - 1L), each = nx)
  a <- id(ci, cj); b <- id(ci + 1L, cj); c <- id(ci + 1L, cj + 1L); d <- id(ci, cj + 1L)
  diag_right <- switch(pattern,
                       right = rep(TRUE, length(a)),
                       left = rep(FALSE, length(a)),
                       alternating = (ci + cj) %% 2L == 0L)
  tris <- rbind(
    cbind(a, b, ifelse(diag_right, c, d)),
    cbind(ifelse(diag_right, a, b), c, d)
  )
  i <- 0:(nx - 1L); j <- 0:(ny - 1L)
  be <- rbind(cbind(id(i, 0L), id(i + 1L, 0L)),
              cbind(id(nx, j), id(nx, j + 1L)),
              cbind(id(i + 1L, ny), id(i, ny)),
              cbind(id(0L, j + 1L), id(0L, j)))
  bm <- c(rep("bottom", nx), rep("right", ny), rep("top", nx), rep("left", ny))
  mesh <- pde_mesh(nodes, tris, be, bm)
  mesh$info$generator <- "rectangle"
  mesh$info$domain <- pde_domain("rectangle", x = x, y = y)
  mesh$info$args <- list(x = x, y = y, nx = nx, ny = ny, pattern = pattern)
  mesh
}

#' Triangular mesh of a polygon
#'
#' Triangulates a simple polygon by ear clipping, improves the coarse
#' triangulation with Lawson edge flips (making it Delaunay while keeping the
#' polygon sides), and then refines it uniformly until the longest edge is at
#' most `h`. Uniform refinement keeps the angles of the coarse
#' triangulation, so mesh quality is set by the polygon's shape.
#'
#' This is a simple, dependency-free mesher for polygonal domains. It is not
#' a quality-guaranteeing constrained Delaunay mesher; for complex geometry,
#' generate the mesh with Gmsh and read it with [read_gmsh()].
#'
#' @param domain A polygon or rectangle [pde_domain()], or a vertex matrix.
#' @param h Target maximum edge length.
#' @param markers Side markers when `domain` is a vertex matrix.
#' @return A `pde_mesh` whose boundary edges carry the polygon side markers.
#' @examples
#' L <- rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
#' mesh <- pde_mesh_polygon(L, h = 0.5)
#' summary(mesh)
#' @export
pde_mesh_polygon <- function(domain, h = Inf, markers = "boundary") {
  if (!inherits(domain, "pde_domain")) {
    domain <- pde_domain("polygon", vertices = domain, markers = markers)
  }
  if (!domain$type %in% c("polygon", "rectangle")) {
    stopf("`domain` must be a polygon or rectangle domain.")
  }
  check_number(h, "h", positive = TRUE, finite = FALSE)
  v <- domain$vertices
  nv <- nrow(v)
  be <- cbind(1:nv, c(2:nv, 1L))
  bm <- domain$markers
  loop <- seq_len(nv)
  if (length(domain$holes)) {
    # Bridge every hole to the boundary so far, giving one weakly simple
    # polygon that ear clipping triangulates.
    hl <- domain$holes
    off <- nv + c(0L, cumsum(vapply(hl, nrow, integer(1))))
    v <- rbind(v, do.call(rbind, hl))
    for (k in seq_along(hl)) {
      ids <- off[k] + seq_len(nrow(hl[[k]]))
      be <- rbind(be, cbind(ids, c(ids[-1], ids[1])))
      bm <- c(bm, rep(domain$hole_markers[k], length(ids)))
    }
    inside <- function(p) points_in_polygon(p, domain$vertices) &&
      !any(vapply(hl, function(h) points_in_polygon(p, h), logical(1)))
    for (k in order(-vapply(hl, function(h) max(h[, 1]), numeric(1)))) {
      ids <- rev(off[k] + seq_len(nrow(hl[[k]])))          # clockwise
      ids <- ids[c(seq(which.max(v[ids, 1]), length(ids)), seq_len(which.max(v[ids, 1]) - 1L))]
      M <- ids[1]
      cand <- unique(loop); cand <- cand[order(rowSums(sweep(v[cand, , drop = FALSE], 2, v[M, ])^2))]
      edges <- be
      P <- NA_integer_
      for (q in cand) {
        ok <- inside((v[q, ] + v[M, ]) / 2)
        if (ok) for (e in seq_len(nrow(edges))) {
          a <- edges[e, 1]; b <- edges[e, 2]
          if (a %in% c(q, M) || b %in% c(q, M)) next
          if (segments_touch(v[q, ], v[M, ], v[a, ], v[b, ])) { ok <- FALSE; break }
        }
        if (ok) { P <- q; break }
      }
      if (is.na(P)) stopf("Could not connect hole %d to the boundary.", k)
      pos <- match(P, loop)
      loop <- c(loop[seq_len(pos)], ids, M, P, loop[-seq_len(pos)])
    }
  }
  tris <- ear_clip(v, loop)
  mesh <- pde_mesh(v, tris, be, bm)
  mesh <- lawson_flip(mesh)
  while (max(mesh_edge_lengths(mesh)) > h) {
    mesh <- refine_uniform(mesh)
  }
  mesh$info$generator <- "polygon"
  mesh$info$domain <- domain
  mesh
}

ear_clip <- function(v, idx = seq_len(nrow(v))) {
  out <- vector("list", length(idx) - 2L)
  k <- 0L
  cross <- function(a, b, c) {
    (v[b, 1] - v[a, 1]) * (v[c, 2] - v[a, 2]) - (v[b, 2] - v[a, 2]) * (v[c, 1] - v[a, 1])
  }
  guard <- 0L
  while (length(idx) > 3L) {
    nn <- length(idx)
    found <- FALSE
    for (p in seq_len(nn)) {
      i0 <- idx[(p - 2L) %% nn + 1L]; i1 <- idx[p]; i2 <- idx[p %% nn + 1L]
      if (cross(i0, i1, i2) <= 0) next
      others <- setdiff(idx, c(i0, i1, i2))
      inside <- FALSE
      for (q in others) {
        if (cross(i0, i1, q) >= 0 && cross(i1, i2, q) >= 0 && cross(i2, i0, q) >= 0) {
          inside <- TRUE; break
        }
      }
      if (inside) next
      k <- k + 1L
      out[[k]] <- c(i0, i1, i2)
      idx <- idx[-p]
      found <- TRUE
      break
    }
    guard <- guard + 1L
    if (!found || guard > 10L * length(out)) stopf("Ear clipping failed; is the polygon simple?")
  }
  k <- k + 1L
  out[[k]] <- idx
  do.call(rbind, out)
}

# Flip interior edges until the triangulation is Delaunay. Boundary edges
# are never flipped, so polygon sides and their markers are preserved.
lawson_flip <- function(mesh, max_pass = 100L) {
  for (pass in seq_len(max_pass)) {
    nodes <- mesh$nodes; tris <- mesh$triangles
    interior <- which(!is.na(mesh$edge_elements[, 2]))
    flipped <- FALSE
    touched <- logical(nrow(tris))
    for (e in interior) {
      t1 <- mesh$edge_elements[e, 1]; t2 <- mesh$edge_elements[e, 2]
      if (touched[t1] || touched[t2]) next
      p <- mesh$edges[e, 1]; q <- mesh$edges[e, 2]
      r1 <- setdiff(tris[t1, ], c(p, q)); r2 <- setdiff(tris[t2, ], c(p, q))
      ang <- function(r) {
        u <- nodes[p, ] - nodes[r, ]; w <- nodes[q, ] - nodes[r, ]
        acos(max(-1, min(1, sum(u * w) / sqrt(sum(u^2) * sum(w^2)))))
      }
      if (ang(r1) + ang(r2) > pi + 1e-12) {
        new1 <- c(r1, r2, p); new2 <- c(r2, r1, q)
        d1 <- tri_det(nodes, matrix(new1, 1)); d2 <- tri_det(nodes, matrix(new2, 1))
        if (d1 <= 0 || d2 <= 0) {
          new1 <- new1[c(1, 3, 2)]; new2 <- new2[c(1, 3, 2)]
          d1 <- -d1; d2 <- -d2
        }
        if (d1 <= 0 || d2 <= 0) next
        tris[t1, ] <- new1; tris[t2, ] <- new2
        touched[c(t1, t2)] <- TRUE
        flipped <- TRUE
      }
    }
    if (!flipped) break
    bd <- mesh$boundary
    m2 <- build_topology(nodes, tris)
    m2$boundary <- assign_boundary(m2, cbind(bd$n1, bd$n2), bd$marker)
    m2$info <- mesh$info
    mesh <- m2
  }
  mesh
}

#' Mesh geometry and quality
#'
#' `mesh_quality()` returns per-element geometric measures: area, longest and
#' shortest edge, minimum and maximum interior angle (degrees), and the
#' normalized shape ratio `2 r / R` (inradius over circumradius, equal to 1
#' for an equilateral triangle and tending to 0 for degenerate ones).
#'
#' @param mesh A [pde_mesh()].
#' @return A data frame with one row per element.
#' @examples
#' q <- mesh_quality(pde_mesh_rectangle(nx = 2))
#' head(q)
#' @export
mesh_quality <- function(mesh) {
  stopifnot(inherits(mesh, "pde_mesh"))
  P <- mesh$nodes; T <- mesh$triangles
  l <- function(a, b) sqrt(rowSums((P[T[, a], , drop = FALSE] - P[T[, b], , drop = FALSE])^2))
  la <- l(2, 3); lb <- l(3, 1); lc <- l(1, 2)
  area <- tri_det(P, T) / 2
  angle <- function(opp, s1, s2) acos(pmin(1, pmax(-1, (s1^2 + s2^2 - opp^2) / (2 * s1 * s2))))
  A1 <- angle(la, lb, lc); A2 <- angle(lb, lc, la); A3 <- pi - A1 - A2
  s <- (la + lb + lc) / 2
  r <- area / s
  R <- la * lb * lc / (4 * area)
  data.frame(
    area = area,
    h_max = pmax(la, lb, lc),
    h_min = pmin(la, lb, lc),
    min_angle = pmin(A1, A2, A3) * 180 / pi,
    max_angle = pmax(A1, A2, A3) * 180 / pi,
    shape_ratio = 2 * r / R
  )
}

mesh_edge_lengths <- function(mesh) {
  P <- mesh$nodes; E <- mesh$edges
  sqrt(rowSums((P[E[, 1], , drop = FALSE] - P[E[, 2], , drop = FALSE])^2))
}

#' Mesh size
#'
#' @param mesh A [pde_mesh()].
#' @return The largest edge length `h` of the mesh.
#' @examples
#' mesh_size(pde_mesh_rectangle(nx = 4))
#' @export
mesh_size <- function(mesh) max(mesh_edge_lengths(mesh))

#' Boundary queries
#'
#' `boundary_markers()` lists the markers present on a mesh.
#' `boundary_nodes()` returns the sorted indices of nodes on boundary edges
#' with the given markers. Selection is by marker only; no floating-point
#' coordinate comparisons are involved.
#'
#' @param mesh A [pde_mesh()].
#' @param marker Character vector of markers. `"boundary"` selects every
#'   boundary edge.
#' @return A character vector of markers, or an integer vector of node indices.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 2)
#' boundary_markers(mesh)
#' boundary_nodes(mesh, "left")
#' @export
boundary_markers <- function(mesh) unique(mesh$boundary$marker)

#' @rdname boundary_markers
#' @export
boundary_nodes <- function(mesh, marker = "boundary") {
  rows <- boundary_rows(mesh, marker)
  sort(unique(c(mesh$boundary$n1[rows], mesh$boundary$n2[rows])))
}

boundary_rows <- function(mesh, marker) {
  if (!is.character(marker) || length(marker) < 1L) {
    stopf("`marker` must be a character vector.")
  }
  known <- unique(mesh$boundary$marker)
  unknown <- setdiff(marker, c(known, "boundary"))
  if (length(unknown)) {
    stopf("Unknown boundary marker(s): %s. Available: %s.",
          paste(unknown, collapse = ", "), paste(c(known, "boundary"), collapse = ", "))
  }
  if ("boundary" %in% marker) return(seq_len(nrow(mesh$boundary)))
  which(mesh$boundary$marker %in% marker)
}

#' Relabel boundary edges
#'
#' Assigns a new marker to the boundary edges whose midpoints satisfy a
#' predicate. This is a convenience for imported meshes without physical
#' markers; selection uses edge midpoints, so it is robust for straight sides
#' as long as `which` is not evaluated exactly on a corner.
#'
#' @param mesh A [pde_mesh()].
#' @param marker New marker name.
#' @param which A function `function(x, y)` returning `TRUE` for edge
#'   midpoints to relabel.
#' @return The mesh with updated markers.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 4)
#' mesh <- set_boundary_marker(mesh, "inlet", function(x, y) x < 1e-8 & y > 0.5)
#' table(mesh$boundary$marker)
#' @export
set_boundary_marker <- function(mesh, marker, which) {
  stopifnot(inherits(mesh, "pde_mesh"))
  check_name(marker, "marker")
  if (!is.function(which)) stopf("`which` must be a function(x, y).")
  P <- mesh$nodes; bd <- mesh$boundary
  mx <- (P[bd$n1, 1] + P[bd$n2, 1]) / 2
  my <- (P[bd$n1, 2] + P[bd$n2, 2]) / 2
  sel <- check_selection(which(mx, my), nrow(bd), "boundary edge")
  mesh$boundary$marker[sel] <- marker
  mesh$info$markers_edited <- TRUE
  mesh
}

#' @export
print.pde_mesh <- function(x, ...) {
  cat("<pde_mesh> 2D triangular mesh\n")
  cat("  nodes:", nrow(x$nodes), " elements:", nrow(x$triangles),
      " edges:", nrow(x$edges), " boundary edges:", nrow(x$boundary), "\n")
  cat("  h =", fmt_num(mesh_size(x)), "\n")
  cat("  markers:", paste(boundary_markers(x), collapse = ", "), "\n")
  r <- unique(mesh_regions(x))
  if (length(r) > 1L) cat("  regions:", paste(r, collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.pde_mesh <- function(object, ...) {
  q <- mesh_quality(object)
  structure(list(
    nodes = nrow(object$nodes), elements = nrow(object$triangles),
    edges = nrow(object$edges), boundary_edges = nrow(object$boundary),
    h_max = max(q$h_max), h_min = min(q$h_min), area = sum(q$area),
    min_angle = min(q$min_angle), max_angle = max(q$max_angle),
    shape_ratio = stats::quantile(q$shape_ratio, c(0, 0.5, 1)),
    markers = table(object$boundary$marker)
  ), class = "summary.pde_mesh")
}

#' @export
print.summary.pde_mesh <- function(x, ...) {
  cat("Triangular mesh summary\n")
  cat("  nodes", x$nodes, " elements", x$elements, " edges", x$edges,
      " boundary edges", x$boundary_edges, "\n")
  cat("  total area", fmt_num(x$area, 8), "\n")
  cat("  edge length: min", fmt_num(x$h_min), " max", fmt_num(x$h_max), "\n")
  cat("  angles (deg): min", fmt_num(x$min_angle), " max", fmt_num(x$max_angle), "\n")
  cat("  shape ratio 2r/R: min", fmt_num(x$shape_ratio[1]), " median",
      fmt_num(x$shape_ratio[2]), " max", fmt_num(x$shape_ratio[3]), "\n")
  cat("  boundary edges by marker:\n")
  for (nm in names(x$markers)) cat("    ", nm, ": ", x$markers[[nm]], "\n", sep = "")
  invisible(x)
}

#' Piecewise coefficients by region
#'
#' `region_coefficient()` builds a coefficient that takes one value per mesh
#' region ([mesh_regions()]), for example the conductivities of different
#' materials. The result is a vectorized `function(x, y)`, so it can be used
#' wherever a coefficient is accepted (diffusion, reaction, source, Young's
#' modulus, error norms). On the mesh it was built from, assembly reads the
#' element values directly; on other meshes (refined or not) the value at a
#' point is that of the region containing it in the original mesh, so the
#' coefficient keeps its geometry under refinement.
#'
#' At points on an interface between regions the value of one adjacent region
#' is returned; quadrature points and the error estimator never evaluate
#' exactly on element edges.
#'
#' @param mesh The [pde_mesh()] whose regions define the pieces.
#' @param values A named numeric vector or list with one value per region
#'   (all regions of the mesh must be covered).
#' @return A function of class `pde_region_coefficient`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' mesh <- set_region(mesh, "right", function(x, y) x > 0.5)
#' k <- region_coefficient(mesh, c(domain = 1, right = 100))
#' prob <- poisson_problem(mesh, diffusion = k, source = 1, boundary = dirichlet_bc(0))
#' range(solve_pde(prob)$u$values)
#' @export
region_coefficient <- function(mesh, values) {
  stopifnot(inherits(mesh, "pde_mesh"))
  values <- unlist(values)
  r <- mesh_regions(mesh)
  if (!is.numeric(values) || is.null(names(values)) || any(!is.finite(values))) {
    stopf("`values` must be finite numbers named by region.")
  }
  miss <- setdiff(unique(r), names(values))
  if (length(miss)) stopf("No value for region(s): %s.", paste(miss, collapse = ", "))
  unknown <- setdiff(names(values), unique(r))
  if (length(unknown)) stopf("Unknown region(s): %s.", paste(unknown, collapse = ", "))
  elem <- as.numeric(values[r])
  f <- function(x, y) {
    e <- locate_points(mesh, x, y)$element
    if (anyNA(e)) stopf("A region coefficient was evaluated outside its mesh.")
    elem[e]
  }
  structure(f, class = c("pde_region_coefficient", "function"), mesh = mesh, element_values = elem,
            values = values)
}

#' @export
print.pde_region_coefficient <- function(x, ...) {
  v <- attr(x, "values")
  cat("<pde_region_coefficient>", paste(names(v), "=", fmt_num(v), collapse = ", "), "\n")
  invisible(x)
}

# Connected components of a mesh (triangles sharing a vertex are connected):
# list(n = number of components, node = component of each node).
mesh_components <- function(mesh) {
  T <- mesh$triangles
  lab <- seq_len(nrow(mesh$nodes))
  repeat {
    m <- pmin(lab[T[, 1]], lab[T[, 2]], lab[T[, 3]])
    new <- lab
    for (k in 1:3) {
      o <- order(m, decreasing = TRUE)  # smallest label wins on repeated nodes
      new[T[o, k]] <- pmin(new[T[o, k]], m[o])
    }
    new <- new[new]  # pointer jumping
    if (identical(new, lab)) break
    lab <- new
  }
  u <- unique(lab)
  list(n = length(u), node = match(lab, u))
}

# Component of each degree of freedom of a space (P2 edge dofs follow the
# element they belong to).
dof_components <- function(space, cc = mesh_components(space$mesh)) {
  comp <- integer(space$ndof)
  ec <- cc$node[space$mesh$triangles[, 1]]
  comp[as.vector(space$cell_dofs)] <- rep(ec, ncol(space$cell_dofs))
  comp
}

# For a mesh `fine` whose triangles each lie inside one triangle of `coarse`
# (for example after refinement), the containing triangle of each fine
# triangle; NULL when the meshes are not nested this way.
nested_parent <- function(fine, coarse) {
  if (identical(fine$nodes, coarse$nodes) && identical(fine$triangles, coarse$triangles)) {
    return(seq_len(nrow(fine$triangles)))
  }
  P <- fine$nodes; T <- fine$triangles
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
  par <- locate_points(coarse, cx, cy)$element
  if (anyNA(par)) return(NULL)
  Q <- coarse$nodes; C <- coarse$triangles
  for (k in 1:3) {
    x1 <- Q[C[par, 1], 1]; y1 <- Q[C[par, 1], 2]
    b11 <- Q[C[par, 2], 1] - x1; b12 <- Q[C[par, 3], 1] - x1
    b21 <- Q[C[par, 2], 2] - y1; b22 <- Q[C[par, 3], 2] - y1
    det <- b11 * b22 - b12 * b21
    dx <- P[T[, k], 1] - x1; dy <- P[T[, k], 2] - y1
    xi <- (b22 * dx - b12 * dy) / det; eta <- (-b21 * dx + b11 * dy) / det
    tol <- 1e-10
    if (any(xi < -tol | eta < -tol | xi + eta > 1 + tol)) return(NULL)
  }
  par
}
