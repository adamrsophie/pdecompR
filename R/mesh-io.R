#' Read a Gmsh mesh
#'
#' Reads two-dimensional triangular meshes in the ASCII Gmsh formats 2.2 and
#' 4.1. Triangles (element type 2) become mesh elements; line elements
#' (type 1) on the boundary carry boundary markers taken from their physical
#' group names (or numbers when names are absent). Other element types
#' (points, quadrilaterals, higher-order elements) are ignored with a
#' message. Boundary edges not covered by a physical line receive the marker
#' `"boundary"`. Nodes not used by triangles are dropped and the remaining
#' nodes renumbered; clockwise triangles are reoriented.
#'
#' @param file Path to a `.msh` file.
#' @return A [pde_mesh()].
#' @examples
#' f <- system.file("extdata", "square.msh", package = "pdecompR")
#' mesh <- read_gmsh(f)
#' table(mesh$boundary$marker)
#' @export
read_gmsh <- function(file) {
  lines <- trimws(readLines(file, warn = FALSE))
  sec <- function(name) {
    s <- which(lines == paste0("$", name)); e <- which(lines == paste0("$End", name))
    if (!length(s)) return(NULL)
    lines[(s[1] + 1L):(e[1] - 1L)]
  }
  fmt <- sec("MeshFormat")
  if (is.null(fmt)) stopf("Not a Gmsh file: no $MeshFormat section.")
  hdr <- strsplit(fmt[1], "\\s+")[[1]]
  version <- as.numeric(hdr[1])
  if (hdr[2] != "0") stopf("Binary Gmsh files are not supported; save as ASCII.")
  pn <- sec("PhysicalNames")
  phys <- character()
  if (!is.null(pn) && length(pn) > 1L) {
    for (l in pn[-1]) {
      p <- regmatches(l, regexec('^(\\d+)\\s+(\\d+)\\s+"(.*)"$', l))[[1]]
      if (length(p) == 4L && p[2] == "1") phys[p[3]] <- p[4]
    }
  }
  num <- function(l) as.numeric(strsplit(l, "\\s+")[[1]])
  if (version < 3) {
    nd <- sec("Nodes")
    nn <- as.integer(nd[1])
    M <- do.call(rbind, lapply(nd[1L + seq_len(nn)], num))
    tags <- M[, 1]; xyz <- M[, 2:4, drop = FALSE]
    el <- sec("Elements")
    ne <- as.integer(el[1])
    tri <- list(); lin <- list(); lphys <- integer(); other <- 0L
    for (l in el[1L + seq_len(ne)]) {
      v <- num(l); type <- v[2]; ntag <- v[3]
      conn <- v[(4 + ntag):length(v)]
      if (type == 2) tri[[length(tri) + 1L]] <- conn[1:3]
      else if (type == 1) { lin[[length(lin) + 1L]] <- conn[1:2]; lphys <- c(lphys, if (ntag > 0) v[4] else 0) }
      else other <- other + 1L
    }
  } else {
    nd <- sec("Nodes")
    h <- num(nd[1]); nblocks <- h[1]
    pos <- 2L
    tags <- numeric(); xyz <- matrix(numeric(), 0, 3)
    for (b in seq_len(nblocks)) {
      bh <- num(nd[pos]); nb <- bh[4]; pos <- pos + 1L
      if (nb > 0) {
        tags <- c(tags, as.numeric(nd[pos:(pos + nb - 1L)]))
        xyz <- rbind(xyz, do.call(rbind, lapply(nd[(pos + nb):(pos + 2L * nb - 1L)], num)))
      }
      pos <- pos + 2L * nb
    }
    ent_phys <- gmsh4_curve_physicals(sec("Entities"))
    el <- sec("Elements")
    h <- num(el[1]); nblocks <- h[1]
    pos <- 2L
    tri <- list(); lin <- list(); lphys <- integer(); other <- 0L
    for (b in seq_len(nblocks)) {
      bh <- num(el[pos]); edim <- bh[1]; etag <- bh[2]; type <- bh[3]; nb <- bh[4]; pos <- pos + 1L
      rows <- if (nb > 0) el[pos:(pos + nb - 1L)] else character()
      pos <- pos + nb
      for (l in rows) {
        v <- num(l)
        if (type == 2) tri[[length(tri) + 1L]] <- v[2:4]
        else if (type == 1) {
          lin[[length(lin) + 1L]] <- v[2:3]
          pt <- ent_phys[as.character(etag)]
          lphys <- c(lphys, if (is.na(pt)) 0 else pt)
        } else other <- other + 1L
      }
    }
  }
  if (!length(tri)) stopf("The file contains no triangles.")
  if (other > 0L) message(sprintf("Ignored %d element(s) that are not triangles or lines.", other))
  if (diff(range(xyz[, 3])) > 0) stopf("The mesh is not planar (z varies); only 2D meshes are supported.")
  T <- do.call(rbind, tri)
  id <- match(T, tags)
  used <- sort(unique(id))
  T2 <- matrix(match(id, used), ncol = 3)
  nodes <- xyz[used, 1:2, drop = FALSE]
  mesh0 <- pde_mesh(nodes, T2, fix_orientation = TRUE)
  bd <- mesh0$boundary
  markers <- rep("boundary", nrow(bd))
  if (length(lin)) {
    L <- matrix(match(match(do.call(rbind, lin), tags), used), ncol = 2)
    n <- nrow(nodes)
    kb <- edge_key(bd$n1, bd$n2, n)
    kl <- edge_key(L[, 1], L[, 2], n)
    pos <- match(kl, kb)
    lab <- ifelse(as.character(lphys) %in% names(phys), phys[as.character(lphys)], as.character(lphys))
    ok <- !is.na(pos) & lphys != 0
    markers[pos[ok]] <- lab[ok]
  }
  mesh <- pde_mesh(nodes, mesh0$triangles, cbind(bd$n1, bd$n2), markers)
  mesh$info$generator <- "gmsh"
  mesh$info$reoriented <- mesh0$info$reoriented
  mesh
}

gmsh4_curve_physicals <- function(ent) {
  out <- numeric()
  if (is.null(ent)) return(out)
  h <- as.numeric(strsplit(ent[1], "\\s+")[[1]])
  np <- h[1]; nc <- h[2]
  rows <- ent[1L + np + seq_len(nc)]
  for (l in rows) {
    v <- as.numeric(strsplit(l, "\\s+")[[1]])
    # tag, minX, minY, minZ, maxX, maxY, maxZ, numPhysicalTags, tags...
    nph <- v[8]
    out[as.character(v[1])] <- if (nph > 0) v[9] else 0
  }
  out
}

#' Write a VTK file
#'
#' Writes a mesh, finite-element function or solution in the legacy ASCII VTK
#' unstructured-grid format, readable by ParaView and VisIt. Point data are
#' the vertex values; for P2 functions only the vertex values are written
#' (the triangles are linear cells).
#'
#' @param x A [pde_mesh()], [fe_function()] or `pde_solution`.
#' @param file Output path.
#' @param name Name of the point data array.
#' @return The file path, invisibly.
#' @examples
#' f <- tempfile(fileext = ".vtk")
#' write_vtk(pde_mesh_rectangle(nx = 2), f)
#' @export
write_vtk <- function(x, file, name = "u") {
  vals <- NULL
  if (inherits(x, "pde_solution")) x <- x$u
  if (inherits(x, "fe_function")) {
    mesh <- x$space$mesh; vals <- x$values[seq_len(nrow(mesh$nodes))]
  } else if (inherits(x, "pde_mesh")) {
    mesh <- x
  } else stopf("`x` must be a mesh, finite-element function or solution.")
  P <- mesh$nodes; T <- mesh$triangles
  out <- c("# vtk DataFile Version 3.0", "pdecompR output", "ASCII", "DATASET UNSTRUCTURED_GRID",
           sprintf("POINTS %d double", nrow(P)),
           sprintf("%.17g %.17g 0", P[, 1], P[, 2]),
           sprintf("CELLS %d %d", nrow(T), 4L * nrow(T)),
           sprintf("3 %d %d %d", T[, 1] - 1L, T[, 2] - 1L, T[, 3] - 1L),
           sprintf("CELL_TYPES %d", nrow(T)), rep("5", nrow(T)))
  if (!is.null(vals)) {
    out <- c(out, sprintf("POINT_DATA %d", nrow(P)), sprintf("SCALARS %s double 1", name),
             "LOOKUP_TABLE default", sprintf("%.17g", vals))
  }
  writeLines(out, file)
  invisible(file)
}
