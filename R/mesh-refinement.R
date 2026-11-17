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
#' [transfer_function()].
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
  m12 <- n + EE[, 1]; m23 <- n + EE[, 2]; m31 <- n + EE[, 3]
  tris <- rbind(cbind(T[, 1], m12, m31), cbind(m12, T[, 2], m23),
                cbind(m31, m23, T[, 3]), cbind(m12, m23, m31))
  bd <- mesh$boundary
  bm <- n + bd$edge
  be <- rbind(cbind(bd$n1, bm), cbind(bm, bd$n2))
  out <- pde_mesh(nodes, tris, be, c(bd$marker, bd$marker), check_overlap = FALSE)
  out <- inherit_regions(out, mesh)
  out$info <- refined_info(mesh, E, nvb = FALSE)
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
  be <- rbind(cbind(bd$n1, bd$n2)[!split, , drop = FALSE],
              cbind(bd$n1, bmid)[split, , drop = FALSE],
              cbind(bmid, bd$n2)[split, , drop = FALSE])
  bmk <- c(bd$marker[!split], bd$marker[split], bd$marker[split])
  out <- pde_mesh(nodes, tris, be, bmk, check_overlap = FALSE)
  out <- inherit_regions(out, mesh)
  out$info <- refined_info(mesh, E[emark, , drop = FALSE], nvb = TRUE)
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
  if (max(abs(mid - fine$nodes[-seq_len(n0), , drop = FALSE]), 0) > 1e-12 * max(1, abs(P0))) return(NULL)
  E
}

refined_info <- function(mesh, parent_edges, nvb) {
  info <- mesh$info
  info$generator <- "refined"
  info$nvb <- nvb
  info$reoriented <- 0L
  info$new_nodes <- unname(parent_edges)
  info$parent_nodes <- nrow(mesh$nodes)
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
  out
}
