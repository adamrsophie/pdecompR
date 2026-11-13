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
