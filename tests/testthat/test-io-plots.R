test_that("Gmsh files are read with physical markers", {
  f <- system.file("extdata", "square.msh", package = "pdecompR")
  expect_message(m <- read_gmsh(f), "Ignored 1")
  expect_equal(nrow(m$triangles), 8)
  expect_setequal(boundary_markers(m), c("bottom", "right", "top", "left"))
  s <- solve_pde(poisson_problem(m, source = 0, boundary = dirichlet_bc(function(x, y) x + y)))
  expect_lt(max(abs(s$u$values - rowSums(m$nodes))), 1e-12)
})

test_that("Gmsh 4.1 files are read", {
  f <- tempfile(fileext = ".msh")
  writeLines(c("$MeshFormat", "4.1 0 8", "$EndMeshFormat",
               "$PhysicalNames", "1", "1 7 \"wall\"", "$EndPhysicalNames",
               "$Entities", "0 1 1 0",
               "1 0 0 0 1 0 0 1 7 0",
               "1 0 0 0 1 1 0 0 0",
               "$EndEntities",
               "$Nodes", "1 4 1 4", "2 1 0 4", "1", "2", "3", "4",
               "0 0 0", "1 0 0", "1 1 0", "0 1 0", "$EndNodes",
               "$Elements", "2 3 1 3",
               "1 1 1 1", "1 1 2",
               "2 1 2 2", "2 1 2 3", "3 1 3 4",
               "$EndElements"), f)
  m <- read_gmsh(f)
  expect_equal(nrow(m$triangles), 2)
  expect_equal(sum(m$boundary$marker == "wall"), 1)
  expect_equal(sum(m$boundary$marker == "boundary"), 3)
})

test_that("VTK output is written", {
  f <- tempfile(fileext = ".vtk")
  V <- fe_space(pde_mesh_rectangle(nx = 2))
  write_vtk(interpolate(function(x, y) x, V), f)
  txt <- readLines(f)
  expect_true(any(grepl("POINTS 9 double", txt)))
  expect_true(any(grepl("CELLS 8 32", txt)))
  expect_true(any(grepl("POINT_DATA 9", txt)))
  expect_error(write_vtk(1, f), "must be")
})

test_that("plots run without errors", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  m <- pde_mesh_rectangle(nx = 4)
  expect_silent(plot(m))
  s <- solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(0)))
  expect_silent(plot(s))
  expect_silent(plot(s$u, type = "contour"))
  expect_silent(plot(s$u, mesh = TRUE))
  expect_silent(plot_indicators(m, estimate_error(s)))
  cs <- convergence_study(poisson_problem(m, source = f_sin, boundary = dirichlet_bc(0)), u_sin,
                          mesh_sizes = c(4, 8))
  expect_silent(plot(cs))
  tr <- solve_pde(heat_problem(m, initial = u_sin, bcs = dirichlet_bc(0), time = c(0, 0.1)), dt = 0.05)
  expect_silent(plot(tr, time = 0.05))
  ad <- adaptive_solve(poisson_problem(m, source = 1, boundary = dirichlet_bc(0)), max_iterations = 2)
  expect_silent(plot(ad))
})

test_that("Gmsh 2.2 elements without tags and with other types are parsed", {
  f <- tempfile(fileext = ".msh")
  writeLines(c("$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$Nodes", "4",
               "1 0 0 0", "2 1 0 0", "3 1 1 0", "4 0 1 0", "$EndNodes", "$Elements", "4",
               "1 15 1 9 1", "2 2 0 1 2 3", "3 2 3 7 7 1 1 3 4", "4 1 0 1 2", "$EndElements"), f)
  expect_message(read_gmsh(f), "Ignored 1")
  m <- suppressMessages(read_gmsh(f))
  expect_equal(nrow(m$triangles), 2L)
  expect_equal(sort(unique(mesh_regions(m))), c("7", "domain"))
  expect_true(all(m$boundary$marker == "boundary"))
})

test_that("VTU field names and PVD file names are escaped and read back", {
  skip_if_not_installed("xml2")
  V <- fe_space(pde_mesh_rectangle(nx = 2))
  u <- interpolate(function(x, y) x + 2 * y, V)
  d <- withr::local_tempdir()
  f <- file.path(d, "special.vtu")
  nm <- "u & v \"test\" <1>"
  write_vtk(u, f, name = nm, format = "ascii")
  doc <- xml2::read_xml(f)
  arr <- xml2::xml_find_first(doc, sprintf("//PointData/DataArray[@Name=%s]", "'u & v \"test\" <1>'"))
  expect_identical(xml2::xml_attr(arr, "Name"), nm)
  vals <- as.numeric(strsplit(trimws(xml2::xml_text(arr)), "\\s+")[[1]])
  expect_equal(vals, u$values, tolerance = 1e-14)
  h <- solve_pde(heat_problem(V$mesh, initial = 0, time = c(0, 0.1)), dt = 0.05)
  p <- file.path(d, "heat & test.pvd")
  write_pvd(h, p, format = "ascii")
  files <- xml2::xml_attr(xml2::xml_find_all(xml2::read_xml(p), "//DataSet"), "file")
  expect_true(all(startsWith(files, "heat & test")))
  expect_true(all(file.exists(file.path(d, files))))
  b <- file.path(d, "binary.vtu")
  write_vtk(u, b, name = nm, format = "binary")
  head <- readLines(b, n = 40, warn = FALSE)
  head <- head[seq_len(grep("<AppendedData", head)[1] - 1L)]
  expect_true(any(grepl("u &amp; v &quot;test&quot; &lt;1&gt;", head, fixed = TRUE)))
})
