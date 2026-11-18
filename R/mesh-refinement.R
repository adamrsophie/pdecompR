#' Mesh refinement
#'
#' `refine_uniform()` splits every triangle into four by joining its edge
#' midpoints (red refinement). Angles are preserved, so shape regularity is
#' inherited from the parent mesh.
#'
#' `refine_mesh()` refines the marked elements by newest vertex bisection and
#' bisects neighbors as needed to keep the mesh conforming (no hanging nodes).
#' Every edge of a marked element is bisected, so marked elements are split
#' into four children. On the first call the longest edge of each triangle is
#' chosen as its refinement edge; afterwards the newest vertex rule applies.
#' Newest vertex bisection produces at most four similarity classes of
#' triangles per initial triangle, so the minimum angle stays bounded away
#' from zero.
#'
#' Both functions preserve boundary markers: a bisected boundary edge yields
#' two edges with the parent's marker. The returned mesh records the parent
#' edge of each new node in `info$new_nodes` (a two-column matrix of node
#' indices), which makes exact transfer of P1 functions possible, see
#' [transfer_function()]. On boundary parts declared curved with
#' [curved_boundary()], new boundary nodes are projected onto the curve
#' instead of staying at the chord midpoint.
#'
#' @param mesh A [pde_mesh()].
#' @param marked Logical vector over elements, or integer element indices.
#' @return The refined `pde_mesh`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 2)
#' nrow(refine_uniform(mesh)$triangles)
#' fine <- refine_mesh(mesh, marked = 1)
#' nrow(fine$triangles)
#' @export
refine_uniform <- function(mesh) {
  stopifnot(inherits(mesh, "pde_mesh"))
  n <- nrow(mesh$nodes); T <- mesh$triangles; EE <- mesh$element_edges
  E <- mesh$edges
  mid <- (mesh$nodes[E[, 1], , drop = FALSE] + mesh$nodes[E[, 2], , drop = FALSE]) / 2
  nodes <- rbind(mesh$nodes, mid)
  pr <- project_curved(nodes, mesh, n + mesh$boundary$edge)
  nodes <- pr$nodes
  m12 <- n + EE[, 1]; m23 <- n + EE[, 2]; m31 <- n + EE[, 3]
  tris <- rbind(cbind(T[, 1], m12, m31), cbind(m12, T[, 2], m23),
                cbind(m31, m23, T[, 3]), cbind(m12, m23, m31))
  bd <- mesh$boundary
  bm <- n + bd$edge
  be <- rbind(cbind(bd$n1, bm), cbind(bm, bd$n2))
  out <- curved_mesh(nodes, tris, be, c(bd$marker, bd$marker), pr$moved)
  out <- carry_curves(out, mesh, be, c(bd$curve, bd$curve))
  parent <- rep(seq_len(nrow(T)), 4L)
  out$regions <- mesh_regions(mesh)[parent]
  out$info <- refined_info(mesh, E, nvb = FALSE, parent)
  out$info$projected_nodes <- pr$moved
  out
}

#' @rdname refine_uniform
#' @references Mitchell, W. F. (1991). Adaptive refinement for arbitrary
#'   finite-element spaces with hierarchical bases. *Journal of Computational
#'   and Applied Mathematics*, 36, 65-78.
#'
#'   Binev, P., Dahmen, W. and DeVore, R. (2004). Adaptive finite element
#'   methods with convergence rates. *Numerische Mathematik*, 97, 219-268.
#' @export
refine_mesh <- function(mesh, marked) {
  stopifnot(inherits(mesh, "pde_mesh"))
  m <- nrow(mesh$triangles)
  marked <- check_indices(marked, m, "marked")
  if (!length(marked)) return(mesh)
  if (!isTRUE(mesh$info$nvb)) mesh <- nvb_initial_labeling(mesh)
  T <- mesh$triangles; EE <- mesh$element_edges; E <- mesh$edges
  n <- nrow(mesh$nodes)
  emark <- logical(nrow(E))
  emark[EE[marked, ]] <- TRUE
  # Closure: an element with any marked edge must bisect its refinement edge.
  repeat {
    need <- (emark[EE[, 2]] | emark[EE[, 3]]) & !emark[EE[, 1]]
    if (!any(need)) break
    emark[EE[need, 1]] <- TRUE
  }
  new_id <- integer(nrow(E))
  new_id[emark] <- n + seq_len(sum(emark))
  mid <- (mesh$nodes[E[emark, 1], , drop = FALSE] + mesh$nodes[E[emark, 2], , drop = FALSE]) / 2
  nodes <- rbind(mesh$nodes, mid)
  r <- emark[EE[, 1]]; s2 <- emark[EE[, 2]]; s3 <- emark[EE[, 3]]
  z1 <- T[, 1]; z2 <- T[, 2]; z3 <- T[, 3]
  mm <- new_id[EE[, 1]]; m23 <- new_id[EE[, 2]]; m31 <- new_id[EE[, 3]]
  keep <- T[!r, , drop = FALSE]
  # Child A = (z3, z1, m), child B = (z2, z3, m).
  A_plain <- r & !s3; A_split <- r & s3
  B_plain <- r & !s2; B_split <- r & s2
  tris <- rbind(
    keep,
    cbind(z3, z1, mm)[A_plain, , drop = FALSE],
    cbind(mm, z3, m31)[A_split, , drop = FALSE],
    cbind(z1, mm, m31)[A_split, , drop = FALSE],
    cbind(z2, z3, mm)[B_plain, , drop = FALSE],
    cbind(mm, z2, m23)[B_split, , drop = FALSE],
    cbind(z3, mm, m23)[B_split, , drop = FALSE]
  )
  bd <- mesh$boundary
  split <- emark[bd$edge]
  bmid <- new_id[bd$edge]
  pr <- project_curved(nodes, mesh, bmid)
  nodes <- pr$nodes
  be <- rbind(cbind(bd$n1, bd$n2)[!split, , drop = FALSE],
              cbind(bd$n1, bmid)[split, , drop = FALSE],
              cbind(bmid, bd$n2)[split, , drop = FALSE])
  bmk <- c(bd$marker[!split], bd$marker[split], bd$marker[split])
  out <- curved_mesh(nodes, tris, be, bmk, pr$moved)
  cur <- bd$curve
  if (!is.null(cur)) out <- carry_curves(out, mesh, be, c(cur[!split], cur[split], cur[split]))
  parent <- c(which(!r), which(A_plain), which(A_split), which(A_split),
              which(B_plain), which(B_split), which(B_split))
  out$regions <- mesh_regions(mesh)[parent]
  out$info <- refined_info(mesh, E[emark, , drop = FALSE], nvb = TRUE, parent)
  out$info$projected_nodes <- pr$moved
  out
}

# Parent edges of the new nodes when `fine` is one refinement of `coarse`
# (same leading nodes, new nodes at midpoints of edges of `coarse`), else
# NULL.
refinement_edges <- function(coarse, fine) {
  inf <- fine$info; P0 <- coarse$nodes; n0 <- nrow(P0)
  E <- inf$new_nodes
  if (is.null(E) || !identical(inf$parent_nodes, n0) || nrow(fine$nodes) != n0 + nrow(E)) return(NULL)
  if (max(abs(fine$nodes[seq_len(n0), , drop = FALSE] - P0), 0) != 0) return(NULL)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  if (!all(key(E[, 1], E[, 2]) %in% key(coarse$edges[, 1], coarse$edges[, 2]))) return(NULL)
  mid <- (P0[E[, 1], , drop = FALSE] + P0[E[, 2], , drop = FALSE]) / 2
  # Nodes projected onto curved boundaries are exempt from the midpoint test.
  chk <- setdiff(seq_len(nrow(E)), inf$projected_nodes - n0)
  if (max(abs(mid[chk, , drop = FALSE] - fine$nodes[n0 + chk, , drop = FALSE]), 0) > 1e-12 * max(abs(P0))) {
    return(NULL)
  }
  E
}

refined_info <- function(mesh, parent_edges, nvb, parent_elements) {
  info <- mesh$info
  info$generator <- "refined"
  info$nvb <- nvb
  info$reoriented <- 0L
  info$new_nodes <- unname(parent_edges)
  info$parent_nodes <- nrow(mesh$nodes)
  # Element of the coarse mesh containing each fine element.
  info$parent_elements <- parent_elements
  info$args <- NULL
  info
}

# Rotate each triangle so that its longest edge is local edge 1 (v1, v2).
nvb_initial_labeling <- function(mesh) {
  P <- mesh$nodes; T <- mesh$triangles
  l <- function(a, b) rowSums((P[T[, a], , drop = FALSE] - P[T[, b], , drop = FALSE])^2)
  L <- cbind(l(1, 2), l(2, 3), l(3, 1))
  k <- max.col(L, ties.method = "first")
  T2 <- T
  s2 <- k == 2L; s3 <- k == 3L
  T2[s2, ] <- T[s2, c(2, 3, 1), drop = FALSE]
  T2[s3, ] <- T[s3, c(3, 1, 2), drop = FALSE]
  bd <- mesh$boundary
  out <- pde_mesh(P, T2, cbind(bd$n1, bd$n2), bd$marker, regions = mesh_regions(mesh))
  out$info <- mesh$info
  out$info$nvb <- TRUE
  out <- carry_curves(out, mesh, cbind(bd$n1, bd$n2), bd$curve)
  out
}

# Refined mesh from midpoint subdivision. Straight subdivision keeps the
# embedding, so the overlap check is skipped; when nodes were moved onto
# curves the projection can fold or overlap elements, and the mesh is
# validated before it is returned.
curved_mesh <- function(nodes, tris, be, markers, moved) {
  if (!length(moved)) return(pde_mesh(nodes, tris, be, markers, check_overlap = FALSE))
  tryCatch(pde_mesh(nodes, tris, be, markers, check_overlap = TRUE), error = function(e) {
    stopf(paste0("Placing the new boundary nodes on their curves gives an invalid mesh (%s). ",
                 "Refine the coarse mesh near the curved boundary first, or relabel those edges with set_boundary_marker(), which drops the curve."),
          conditionMessage(e))
  })
}

# Curves are attached to boundary edges (`boundary$curve`, an index into
# `mesh$curves`). Copy them to the boundary edges `be` of a mesh derived from
# `mesh`, where `cur` holds the curve of each row of `be`.
carry_curves <- function(out, mesh, be, cur) {
  if (is.null(mesh$curves) || is.null(cur)) return(out)
  n <- nrow(out$nodes)
  pos <- match(edge_key(be[, 1], be[, 2], n), edge_key(out$boundary$n1, out$boundary$n2, n))
  out$boundary$curve <- NA_integer_
  out$boundary$curve[pos] <- cur
  out$curves <- mesh$curves
  out
}

# Move the new nodes `ids` (one per boundary edge of `mesh`, NA or 0 when
# the edge is not split) of curved boundary edges onto their curve.
project_curved <- function(nodes, mesh, ids) {
  cur <- mesh$boundary$curve
  moved <- integer(0)
  if (is.null(cur)) return(list(nodes = nodes, moved = moved))
  for (cv in sort(unique(cur[!is.na(cur)]))) {
    k <- ids[which(cur == cv)]
    k <- k[!is.na(k) & k > 0L]
    if (!length(k)) next
    nodes[k, ] <- eval_projection(mesh$curves[[cv]], nodes[k, 1], nodes[k, 2], names(mesh$curves)[cv])
    moved <- c(moved, k)
  }
  list(nodes = nodes, moved = sort(moved))
}

eval_projection <- function(f, x, y, marker) {
  p <- f(x, y)
  if (!is.matrix(p) || !is.numeric(p) || nrow(p) != length(x) || ncol(p) != 2L || any(!is.finite(p))) {
    stopf("The projection of curved boundary \"%s\" must return a finite two-column matrix with one row per point.", marker)
  }
  p
}

#' Curved boundaries
#'
#' Declares the boundary edges with marker `marker` as approximations of a
#' curve, given by its closest point `projection(x, y)`. Later calls of
#' [refine_uniform()] and [refine_mesh()] then place new boundary nodes on
#' the curve instead of at chord midpoints, so the polygonal domain
#' converges to the curved one and P1 elements keep their second-order L2
#' accuracy (a mesh that keeps the coarse chords stalls at the geometric
#' error). The elements stay straight-sided: P2 edge nodes lie on the
#' chords, which limits P2 to the accuracy of the polygonal domain
#' (isoparametric elements are not provided).
#'
#' The existing nodes on the marked edges must already lie on the curve. The
#' curve is attached to those edges (column `curve` of `mesh$boundary`, an
#' index into `mesh$curves`) and passed to their children by refinement;
#' relabeling edges with [set_boundary_marker()] makes them straight again,
#' whatever marker the curve was declared with. Because refined meshes no longer cover the same
#' domain as their parents, [errornorm()] against a finite-element reference
#' from another mesh of the sequence and [l2_project()] between them stop
#' with an error; compare with exact functions instead. Multigrid
#' prolongation and the one-step P1 path of [transfer_function()] use the
#' chord average at projected nodes; other transfers between meshes of a
#' curved sequence (P2, or across several refinements) stop with an error
#' because the projected nodes lie outside the coarser mesh.
#'
#' @param mesh A [pde_mesh()].
#' @param marker Boundary marker of the curved edges.
#' @param projection `function(x, y)` returning a two-column matrix of the
#'   closest points on the curve.
#' @return The mesh, with the curve stored in `curves` and attached to the
#'   marked edges.
#' @examples
#' n <- 16
#' circle <- cbind(cos(2 * pi * (1:n) / n), sin(2 * pi * (1:n) / n))
#' disk <- curved_boundary(pde_mesh_polygon(circle), "boundary",
#'                         function(x, y) cbind(x, y) / sqrt(x^2 + y^2))
#' fine <- refine_uniform(refine_uniform(disk))
#' range(sqrt(rowSums(fine$nodes[boundary_nodes(fine), ]^2)))
#' @export
curved_boundary <- function(mesh, marker, projection) {
  stopifnot(inherits(mesh, "pde_mesh"))
  check_name(marker, "marker")
  if (!is.function(projection)) stopf("`projection` must be a function(x, y).")
  rows <- boundary_rows(mesh, marker)
  bd <- mesh$boundary
  nd <- sort(unique(c(bd$n1[rows], bd$n2[rows])))
  P <- mesh$nodes[nd, , drop = FALSE]
  q <- eval_projection(projection, P[, 1], P[, 2], marker)
  h <- min(mesh_edge_lengths(mesh))
  if (max(abs(q - P)) > 1e-8 * max(h, 1e-300)) {
    stopf("The nodes on boundary \"%s\" do not lie on the curve (largest distance %s).", marker,
          fmt_num(max(sqrt(rowSums((q - P)^2)))))
  }
  mesh$curves <- c(mesh$curves, stats::setNames(list(projection), marker))
  if (is.null(bd$curve)) mesh$boundary$curve <- NA_integer_
  mesh$boundary$curve[rows] <- length(mesh$curves)
  mesh
}
