#' Read a Gmsh mesh
#'
#' Reads two-dimensional triangular meshes in the ASCII Gmsh formats 2.2 and
#' 4.1. Triangles (element type 2) become mesh elements; line elements
#' (type 1) on the boundary carry boundary markers taken from their physical
#' group names (or numbers when names are absent); physical surfaces become
#' the regions of the triangles ([mesh_regions()]), `"domain"` when absent. Other element types
#' (points, quadrilaterals, higher-order elements) are ignored with a
#' message. Boundary edges not covered by a physical line receive the marker
#' `"boundary"`. Nodes not used by triangles are dropped and the remaining
#' nodes renumbered; clockwise triangles are reoriented.
#'
#' @param file Path to a `.msh` file.
#' @param check_overlap Passed to [pde_mesh()].
#' @return A [pde_mesh()].
#' @examples
#' f <- system.file("extdata", "square.msh", package = "pdecompR")
#' mesh <- read_gmsh(f)
#' table(mesh$boundary$marker)
#' @export
read_gmsh <- function(file, check_overlap = TRUE) {
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
  phys <- character(); phys2 <- character()
  if (!is.null(pn) && length(pn) > 1L) {
    for (l in pn[-1]) {
      p <- regmatches(l, regexec('^(\\d+)\\s+(\\d+)\\s+"(.*)"$', l))[[1]]
      if (length(p) == 4L && p[2] == "1") phys[p[3]] <- p[4]
      if (length(p) == 4L && p[2] == "2") phys2[p[3]] <- p[4]
    }
  }
  num <- function(l) as.numeric(strsplit(l, "\\s+")[[1]])
  # All numbers of many lines at once, with the offset of each line's first
  # number in the flat vector.
  flat <- function(l) {
    sp <- strsplit(l, "\\s+")
    list(v = as.numeric(unlist(sp, use.names = FALSE)), start = cumsum(c(0L, lengths(sp)[-length(sp)])))
  }
  if (version < 3) {
    nd <- sec("Nodes")
    nn <- as.integer(nd[1])
    fn <- flat(nd[1L + seq_len(nn)])
    tags <- fn$v[fn$start + 1]
    xyz <- cbind(fn$v[fn$start + 2], fn$v[fn$start + 3], fn$v[fn$start + 4])
    el <- sec("Elements")
    ne <- as.integer(el[1])
    fe <- flat(el[1L + seq_len(ne)])
    st <- fe$start; type <- fe$v[st + 2]; ntag <- fe$v[st + 3]
    ephys <- ifelse(ntag > 0, fe$v[st + 4], 0)
    c0 <- st + 4 + ntag
    it <- type == 2; il <- type == 1
    tri <- cbind(fe$v[c0[it]], fe$v[c0[it] + 1], fe$v[c0[it] + 2]); tphys <- ephys[it]
    lin <- cbind(fe$v[c0[il]], fe$v[c0[il] + 1]); lphys <- ephys[il]
    other <- sum(!it & !il)
  } else {
    nd <- sec("Nodes")
    h <- num(nd[1]); nblocks <- h[1]
    pos <- 2L
    tag_rows <- integer(); xyz_rows <- integer()
    for (b in seq_len(nblocks)) {
      bh <- num(nd[pos]); nb <- bh[4]; pos <- pos + 1L
      if (bh[3] != 0) stopf("Parametric node blocks in Gmsh 4.1 files are not supported; save without parametric coordinates.")
      if (nb > 0) {
        tag_rows <- c(tag_rows, pos:(pos + nb - 1L))
        xyz_rows <- c(xyz_rows, (pos + nb):(pos + 2L * nb - 1L))
      }
      pos <- pos + 2L * nb
    }
    tags <- as.numeric(nd[tag_rows])
    fx <- flat(nd[xyz_rows])
    xyz <- cbind(fx$v[fx$start + 1], fx$v[fx$start + 2], fx$v[fx$start + 3])
    ent_phys <- gmsh4_curve_physicals(sec("Entities"))
    surf_phys <- gmsh4_curve_physicals(sec("Entities"), dim = 2L)
    el <- sec("Elements")
    h <- num(el[1]); nblocks <- h[1]
    pos <- 2L
    rows <- integer(); btype <- numeric(); bphys <- numeric()
    for (b in seq_len(nblocks)) {
      bh <- num(el[pos]); etag <- bh[2]; type <- bh[3]; nb <- bh[4]; pos <- pos + 1L
      if (nb > 0) {
        rows <- c(rows, pos:(pos + nb - 1L))
        pt <- if (type == 2) surf_phys[as.character(etag)] else if (type == 1) ent_phys[as.character(etag)] else NA
        btype <- c(btype, rep(type, nb)); bphys <- c(bphys, rep(if (is.na(pt)) 0 else pt, nb))
      }
      pos <- pos + nb
    }
    fe <- flat(el[rows]); st <- fe$start
    it <- btype == 2; il <- btype == 1
    tri <- cbind(fe$v[st[it] + 2], fe$v[st[it] + 3], fe$v[st[it] + 4]); tphys <- unname(bphys[it])
    lin <- cbind(fe$v[st[il] + 2], fe$v[st[il] + 3]); lphys <- unname(bphys[il])
    other <- sum(!it & !il)
  }
  if (!nrow(tri)) {
    stopf("The file contains no linear triangles (element type 2); second-order meshes are not supported.")
  }
  if (any(!match(tri, tags, nomatch = 0L))) stopf("A triangle references a node that is not in the file.")
  if (other > 0L) message(sprintf("Ignored %d element(s) that are not triangles or lines.", other))
  if (diff(range(xyz[, 3])) > 0) stopf("The mesh is not planar (z varies); only 2D meshes are supported.")
  T <- tri
  id <- match(T, tags)
  used <- sort(unique(id))
  T2 <- matrix(match(id, used), ncol = 3)
  nodes <- xyz[used, 1:2, drop = FALSE]
  mesh0 <- pde_mesh(nodes, T2, fix_orientation = TRUE, check_overlap = check_overlap)
  bd <- mesh0$boundary
  markers <- rep("boundary", nrow(bd))
  if (nrow(lin)) {
    L <- matrix(match(match(lin, tags), used), ncol = 2)
    n <- nrow(nodes)
    kb <- edge_key(bd$n1, bd$n2, n)
    kl <- edge_key(L[, 1], L[, 2], n)
    pos <- match(kl, kb)
    lab <- ifelse(as.character(lphys) %in% names(phys), phys[as.character(lphys)], as.character(lphys))
    ok <- !is.na(pos) & lphys != 0
    markers[pos[ok]] <- lab[ok]
  }
  regions <- ifelse(as.character(tphys) %in% names(phys2), phys2[as.character(tphys)],
                    ifelse(tphys == 0, "domain", as.character(tphys)))
  mesh <- pde_mesh(nodes, mesh0$triangles, cbind(bd$n1, bd$n2), markers, regions = unname(regions),
                   check_overlap = FALSE)
  mesh$info$generator <- "gmsh"
  mesh$info$reoriented <- mesh0$info$reoriented
  mesh
}

gmsh4_curve_physicals <- function(ent, dim = 1L) {
  out <- numeric()
  if (is.null(ent)) return(out)
  h <- as.numeric(strsplit(ent[1], "\\s+")[[1]])
  np <- h[1]; nc <- h[2]; ns <- h[3]
  rows <- if (dim == 1L) ent[1L + np + seq_len(nc)] else ent[1L + np + nc + seq_len(ns)]
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
#' unstructured-grid format, readable by ParaView and VisIt. P1 data are
#' written on linear triangles; P2 data on quadratic triangles (VTK cell type
#' 22) with all six nodal values, which ParaView renders with curved
#' interpolation. A displacement solution of [elasticity_problem()] is written
#' as a `VECTORS` array, optionally with the nodal stresses of
#' [elasticity_stress()] (on the vertices, as P1 data). A time-dependent
#' solution is written at one stored time (all species for reaction-diffusion
#' systems); [write_pvd()] writes a whole series.
#'
#' A `file` ending in `.vtu` is written in the XML VTK format instead, with
#' the same content: by default as raw little-endian binary appended data
#' (`format = "binary"`; Float64 values, Int32 connectivity, UInt32 block
#' headers), which is exact and faster to write and read than text, or as
#' text with `format = "ascii"`.
#'
#' @param x A [pde_mesh()], [fe_function()], `pde_solution`,
#'   `pde_vector_solution` or `pde_transient`.
#' @param file Output path (`.vtk` for legacy, `.vtu` for XML).
#' @param name Name of the point data array (scalar data and displacements).
#' @param time Stored time to write for a `pde_transient` (default: the last).
#' @param stress For elasticity solutions, also write `sxx`, `syy`, `sxy` and
#'   `von_mises` (P1 solutions only, since the stresses are vertex values).
#' @param format Encoding of `.vtu` files: `"binary"` or `"ascii"`.
#' @return The file path, invisibly.
#' @examples
#' f <- tempfile(fileext = ".vtk")
#' write_vtk(pde_mesh_rectangle(nx = 2), f)
#' @export
write_vtk <- function(x, file, name = "u", time = NULL, stress = FALSE, format = c("binary", "ascii")) {
  format <- match.arg(format)
  d <- vtk_data(x, name, time, stress)
  if (grepl("\\.vtu$", file, ignore.case = TRUE)) {
    if (format == "binary") write_vtu_binary(d, file, name) else write_vtu_file(d, file, name)
  } else write_legacy_vtk(d, file, name)
  invisible(file)
}

vtk_data <- function(x, name, time, stress) {
  arrays <- list()
  vectors <- NULL
  space <- NULL
  if (inherits(x, "pde_transient")) {
    if (is.null(time)) time <- max(x$times)
    sp <- if (is.list(x$values)) names(x$values) else name
    arrays <- lapply(stats::setNames(seq_along(sp), sp), function(i) solution_at(x, time, i)$values)
    space <- x$space
  } else if (inherits(x, "pde_vector_solution")) {
    space <- x$space
    vectors <- cbind(x$ux$values, x$uy$values)
    if (stress) {
      if (space$degree != 1L) stopf("`stress = TRUE` is available for P1 displacements only.")
      st <- elasticity_stress(x, "nodes")
      arrays <- lapply(st[c("sxx", "syy", "sxy", "von_mises")], `[[`, "values")
    }
  } else {
    if (inherits(x, "pde_solution")) x <- x$u
    if (inherits(x, "fe_function")) {
      space <- x$space; arrays[[name]] <- x$values
    } else if (inherits(x, "pde_mesh")) {
      space <- fe_space(x, 1)
    } else stopf("`x` must be a mesh, finite-element function or solution.")
  }
  quad <- space$degree == 2L
  list(P = if (quad) space$dof_coords else space$mesh$nodes,
       C = if (quad) space$cell_dofs else space$mesh$triangles,
       type = if (quad) 22L else 5L, arrays = arrays, vectors = vectors)
}

write_legacy_vtk <- function(d, file, name) {
  P <- d$P; C <- d$C; arrays <- d$arrays; vectors <- d$vectors
  np <- nrow(P); nc <- nrow(C); k <- ncol(C)
  out <- c("# vtk DataFile Version 3.0", "pdecompR output", "ASCII", "DATASET UNSTRUCTURED_GRID",
           sprintf("POINTS %d double", np),
           sprintf("%.17g %.17g 0", P[, 1], P[, 2]),
           sprintf("CELLS %d %d", nc, (k + 1L) * nc),
           do.call(paste, c(list(k), as.data.frame(C - 1L))),
           sprintf("CELL_TYPES %d", nc), rep(as.character(d$type), nc))
  if (length(arrays) || !is.null(vectors)) out <- c(out, sprintf("POINT_DATA %d", np))
  for (nm in names(arrays)) {
    out <- c(out, sprintf("SCALARS %s double 1", nm), "LOOKUP_TABLE default",
             sprintf("%.17g", arrays[[nm]][seq_len(np)]))
  }
  if (!is.null(vectors)) {
    out <- c(out, sprintf("VECTORS %s double", name),
             sprintf("%.17g %.17g 0", vectors[, 1], vectors[, 2]))
  }
  writeLines(out, file)
}

write_vtu_file <- function(d, file, name) {
  P <- d$P; C <- d$C; np <- nrow(P); nc <- nrow(C); k <- ncol(C)
  da <- function(nm, comps, v) c(sprintf('<DataArray type="Float64" Name="%s" NumberOfComponents="%d" format="ascii">', xml_escape(nm), comps),
                                 v, "</DataArray>")
  pd <- character(0)
  for (nm in names(d$arrays)) pd <- c(pd, da(nm, 1L, sprintf("%.17g", d$arrays[[nm]][seq_len(np)])))
  if (!is.null(d$vectors)) pd <- c(pd, da(name, 3L, sprintf("%.17g %.17g 0", d$vectors[, 1], d$vectors[, 2])))
  out <- c('<?xml version="1.0"?>',
           '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
           "<UnstructuredGrid>", sprintf('<Piece NumberOfPoints="%d" NumberOfCells="%d">', np, nc),
           "<PointData>", pd, "</PointData>",
           "<Points>", da("Points", 3L, sprintf("%.17g %.17g 0", P[, 1], P[, 2])), "</Points>",
           "<Cells>",
           '<DataArray type="Int64" Name="connectivity" format="ascii">',
           do.call(paste, as.data.frame(C - 1L)), "</DataArray>",
           '<DataArray type="Int64" Name="offsets" format="ascii">', paste(k * seq_len(nc), collapse = " "),
           "</DataArray>",
           '<DataArray type="UInt8" Name="types" format="ascii">', paste(rep(d$type, nc), collapse = " "),
           "</DataArray>", "</Cells>", "</Piece>", "</UnstructuredGrid>", "</VTKFile>")
  writeLines(out, file)
}

# XML VTK with raw appended data: every array is a UInt32 byte count followed
# by its little-endian values; DataArray offsets point into that block.
write_vtu_binary <- function(d, file, name) {
  P <- d$P; C <- d$C; np <- nrow(P); nc <- nrow(C); k <- ncol(C)
  blocks <- list(); tags <- character(0); off <- 0
  add <- function(nm, type, comps, v) {
    size <- switch(type, Float64 = 8L, Int32 = 4L, UInt8 = 1L)
    con <- rawConnection(raw(0), "wb")
    writeBin(as.integer(length(v) * size), con, size = 4L, endian = "little")
    if (type == "Float64") writeBin(as.double(v), con, size = 8L, endian = "little")
    else writeBin(as.integer(v), con, size = size, endian = "little")
    b <- rawConnectionValue(con); close(con)
    tag <- sprintf('<DataArray type="%s" Name="%s" NumberOfComponents="%d" format="appended" offset="%.0f"/>',
                   type, xml_escape(nm), comps, off)
    blocks[[length(blocks) + 1L]] <<- b; off <<- off + length(b)
    tag
  }
  xyz <- function(M) as.vector(t(cbind(M[, 1], M[, 2], 0)))
  pd <- character(0)
  for (nm in names(d$arrays)) pd <- c(pd, add(nm, "Float64", 1L, d$arrays[[nm]][seq_len(np)]))
  if (!is.null(d$vectors)) pd <- c(pd, add(name, "Float64", 3L, xyz(d$vectors)))
  pts <- add("Points", "Float64", 3L, xyz(P))
  cells <- c(add("connectivity", "Int32", 1L, as.vector(t(C - 1L))), add("offsets", "Int32", 1L, k * seq_len(nc)),
             add("types", "UInt8", 1L, rep(d$type, nc)))
  head <- c('<?xml version="1.0"?>',
            '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">',
            "<UnstructuredGrid>", sprintf('<Piece NumberOfPoints="%d" NumberOfCells="%d">', np, nc),
            "<PointData>", pd, "</PointData>", "<Points>", pts, "</Points>", "<Cells>", cells, "</Cells>",
            "</Piece>", "</UnstructuredGrid>", '<AppendedData encoding="raw">')
  con <- file(file, "wb"); on.exit(close(con))
  writeBin(charToRaw(paste0(paste(head, collapse = "\n"), "\n_")), con)
  for (b in blocks) writeBin(b, con)
  writeBin(charToRaw("\n</AppendedData>\n</VTKFile>\n"), con)
}

#' Write a time series for ParaView
#'
#' Writes every stored time of a time-dependent solution as an XML VTK file
#' (`.vtu`) and a ParaView collection file (`.pvd`) that lists them with
#' their times, so ParaView opens the series as one animated data set.
#'
#' @param x A `pde_transient` from [solve_pde()].
#' @param file Path of the `.pvd` file; the `.vtu` files are written next to
#'   it with the same stem and the time index appended.
#' @param name Name of the point data array for single-species solutions.
#' @param format Encoding of the `.vtu` files, as in [write_vtk()].
#' @return The `.pvd` path, invisibly.
#' @examples
#' p <- heat_problem(pde_mesh_rectangle(nx = 4), initial = function(x, y) x * y,
#'                   bcs = dirichlet_bc(0), time = c(0, 0.1))
#' s <- solve_pde(p, dt = 0.05)
#' f <- tempfile(fileext = ".pvd")
#' write_pvd(s, f)
#' @export
write_pvd <- function(x, file, name = "u", format = c("binary", "ascii")) {
  format <- match.arg(format)
  if (!inherits(x, "pde_transient")) stopf("`x` must be a time-dependent solution.")
  stem <- sub("\\.pvd$", "", file, ignore.case = TRUE)
  files <- sprintf("%s_%04d.vtu", stem, seq_along(x$times))
  for (i in seq_along(x$times)) write_vtk(x, files[i], name = name, time = x$times[i], format = format)
  ds <- sprintf('<DataSet timestep="%.17g" part="0" file="%s"/>', x$times, xml_escape(basename(files)))
  writeLines(c('<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1">', "<Collection>", ds,
               "</Collection>", "</VTKFile>"), file)
  invisible(file)
}

# Escape text for an XML attribute value.
xml_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE); gsub("'", "&apos;", x, fixed = TRUE)
}
