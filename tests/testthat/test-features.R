test_that("direct solvers flag singular matrices and nonlinear pure Neumann problems stop", {
  V <- fe_space(pde_mesh_rectangle(nx = 8))
  K <- assemble_matrix(V, "stiffness")
  b <- assemble_vector(V, function(x, y) x - 0.5)
  for (m in c("auto", "cholesky", "lu")) expect_equal(solve_linear_system(K, b, m)$info$status, "near-singular")
  expect_equal(solve_linear_system(K + assemble_matrix(V, "mass"), b)$info$status, "solved")
  expect_error(solve_pde(nonlinear_problem(pde_mesh_rectangle(nx = 4), diffusion = function(u) 1 + u^2,
                                           source = 0)), "not unique")
})

test_that("elasticity stresses are exact for a uniaxial bar", {
  m <- pde_mesh_rectangle(c(0, 2), c(0, 1), nx = 8, ny = 4)
  bcs <- list(displacement_bc(0, "left", components = "x"), displacement_bc(0, "bottom", components = "y"),
              traction_bc(c(1, 0), "right"))
  s <- solve_pde(elasticity_problem(m, E = 1, nu = 0.25, plane = "stress", bcs = bcs))
  st <- elasticity_stress(s)
  expect_equal(st$sxx, rep(1, nrow(m$triangles)), tolerance = 1e-10)
  expect_lt(max(abs(st$syy), abs(st$sxy), abs(st$szz)), 1e-10)
  expect_equal(st$von_mises, rep(1, nrow(m$triangles)), tolerance = 1e-10)
  expect_equal(st$exx, rep(1, nrow(m$triangles)), tolerance = 1e-10)
  nd <- elasticity_stress(s, "nodes")
  expect_s3_class(nd$von_mises, "fe_function")
  s2 <- solve_pde(elasticity_problem(m, E = 1, nu = 0.25, plane = "strain", bcs = bcs))
  st2 <- elasticity_stress(s2)
  expect_equal(st2$szz, 0.25 * (st2$sxx + st2$syy), tolerance = 1e-10)
  f <- tempfile(fileext = ".vtk")
  write_vtk(s, f, stress = TRUE)
  expect_true(any(grepl("^SCALARS von_mises", readLines(f))))
})

test_that("P2 functions are written as quadratic VTK cells", {
  V <- fe_space(pde_mesh_rectangle(nx = 2), 2)
  f <- tempfile(fileext = ".vtk")
  write_vtk(interpolate(function(x, y) x * y, V), f)
  l <- readLines(f)
  expect_true(any(l == sprintf("POINTS %d double", V$ndof)))
  expect_equal(l[grep("^CELL_TYPES", l) + 1L], "22")
})

test_that("transient solutions record the integral of the solution", {
  m <- pde_mesh_rectangle(nx = 8)
  s <- solve_pde(heat_problem(m, initial = function(x, y) x, time = c(0, 0.5)), dt = 0.1)
  expect_equal(s$diagnostics$integral, rep(0.5, 6), tolerance = 1e-12)
  r <- solve_pde(reaction_diffusion_problem(m, species = c("a", "b"),
                                            reaction = function(u, x, y, t) cbind(-u[, 1] * u[, 2], u[, 1] * u[, 2]),
                                            initial = list(a = 1, b = function(x, y) x)), dt = 0.1)
  expect_equal(rowSums(r$diagnostics$integral), rep(1.5, 11), tolerance = 1e-10)
})

test_that("time-dependent coefficients keep the temporal orders", {
  ex <- function(x, y, t) exp(-t) * (x + 2 * y)
  p <- parabolic_problem(pde_mesh_rectangle(nx = 4), diffusion = function(x, y, t) 1 + t + 0 * x,
                         reaction = function(x, y, t) t + 0 * x,
                         source = function(x, y, t) (t - 1) * ex(x, y, t),
                         initial = function(x, y) x + 2 * y, bcs = dirichlet_bc(ex), time = c(0, 1))
  err <- function(dt, method) {
    s <- solve_pde(p, dt = dt, method = method, output_times = 1)
    xy <- s$space$dof_coords
    max(abs(s$values[, 1] - ex(xy[, 1], xy[, 2], 1)))
  }
  e_cn <- vapply(c(0.05, 0.025), err, numeric(1), method = "crank-nicolson")
  e_be <- vapply(c(0.05, 0.025), err, numeric(1), method = "backward-euler")
  expect_gt(log2(e_cn[1] / e_cn[2]), 1.9)
  expect_gt(log2(e_be[1] / e_be[2]), 0.9)
  expect_error(parabolic_problem(pde_mesh_rectangle(nx = 2), diffusion = function(x, y, t) 1 + t + 0 * x,
                                 advection = c(1, 0), stabilization = "supg") |>
                 solve_pde(dt = 0.5), "SUPG")
})

test_that("regions carry piecewise coefficients through refinement", {
  m <- set_region(pde_mesh_rectangle(nx = 8), "right", function(x, y) x > 0.5)
  expect_equal(sort(unique(mesh_regions(m))), c("domain", "right"))
  k <- region_coefficient(m, c(domain = 1, right = 100))
  expect_output(print(k), "right = +100")
  ex <- function(x, y) ifelse(x < 0.5, x * 200 / 101, 100 / 101 + (x - 0.5) * 2 / 101)
  bcs <- list(dirichlet_bc(0, "left"), dirichlet_bc(1, "right"))
  s <- solve_pde(elliptic_problem(m, diffusion = k, source = 0, bcs = bcs))
  expect_lt(errornorm(s, ex), 1e-12)
  expect_lt(estimate_error(s)$total, 1e-8)
  r <- refine_uniform(m)
  expect_equal(as.numeric(table(mesh_regions(r))), c(256, 256))
  expect_lt(errornorm(solve_pde(elliptic_problem(r, diffusion = k, source = 0, bcs = bcs)), ex), 1e-12)
  loc <- refine_mesh(m, c(1, 70))
  P <- loc$nodes; T <- loc$triangles
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  expect_true(all((cx > 0.5) == (mesh_regions(loc) == "right")))
  expect_error(region_coefficient(m, c(domain = 1)), "No value")
  expect_error(region_coefficient(m, c(domain = 1, right = 2, other = 3)), "Unknown")
})

test_that("Gmsh physical surfaces become regions", {
  f <- tempfile(fileext = ".msh")
  writeLines(c("$MeshFormat", "2.2 0 8", "$EndMeshFormat",
               "$PhysicalNames", "3", "1 1 \"wall\"", "2 2 \"left\"", "2 3 \"right\"", "$EndPhysicalNames",
               "$Nodes", "6", "1 0 0 0", "2 1 0 0", "3 2 0 0", "4 0 1 0", "5 1 1 0", "6 2 1 0", "$EndNodes",
               "$Elements", "5", "1 1 2 1 1 1 2", "2 2 2 2 1 1 2 5", "3 2 2 2 1 1 5 4",
               "4 2 2 3 2 2 3 6", "5 2 2 3 2 2 6 5", "$EndElements"), f)
  g <- read_gmsh(f)
  expect_equal(mesh_regions(g), c("left", "left", "right", "right"))
})

test_that("L2 projection is conservative across unrelated meshes", {
  a <- interpolate(function(x, y) sin(5 * x) * y, fe_space(pde_mesh_rectangle(nx = 8)))
  I <- function(u) sum(assemble_matrix(u$space, "mass") %*% u$values)
  for (V in list(fe_space(pde_mesh_rectangle(nx = 5, pattern = "left")),
                 fe_space(pde_mesh_rectangle(nx = 3), 2), fe_space(refine_uniform(a$space$mesh)))) {
    expect_equal(I(l2_project(a, V)), I(a), tolerance = 1e-13)
  }
  g <- function(x, y) 1 + 2 * x - y
  q <- l2_project(interpolate(g, fe_space(pde_mesh_rectangle(nx = 3))), fe_space(pde_mesh_rectangle(nx = 4, pattern = "left")))
  expect_lt(max(abs(q$values - g(q$space$dof_coords[, 1], q$space$dof_coords[, 2]))), 1e-12)
  expect_equal(I(l2_project(g, fe_space(pde_mesh_rectangle(nx = 4)))), 1.5, tolerance = 1e-12)
})

test_that("VTU and PVD output", {
  p <- heat_problem(pde_mesh_rectangle(nx = 2), initial = function(x, y) x * y, bcs = dirichlet_bc(0), time = c(0, 0.1))
  s <- solve_pde(p, dt = 0.05)
  d <- tempfile(); dir.create(d)
  write_pvd(s, file.path(d, "a.pvd"))
  expect_setequal(list.files(d), c("a.pvd", sprintf("a_%04d.vtu", 1:3)))
  x <- readLines(file.path(d, "a_0001.vtu"))
  expect_true(any(grepl('NumberOfPoints="9" NumberOfCells="8"', x, fixed = TRUE)))
  expect_equal(sum(grepl('timestep=', readLines(file.path(d, "a.pvd")))), 3L)
  f <- file.path(d, "q.vtu")
  write_vtk(fe_function(fe_space(pde_mesh_rectangle(nx = 1), 2), 1:9), f)
  y <- readLines(f)
  expect_equal(y[which(y == '<DataArray type="UInt8" Name="types" format="ascii">') + 1L], "22 22")
  expect_error(write_pvd(pde_mesh_rectangle(nx = 1), f), "time-dependent")
})

test_that("spatially varying diffusion tensors", {
  m <- pde_mesh_rectangle(nx = 8)
  a <- solve_pde(elliptic_problem(m, diffusion = matrix(c(2, 0.3, 0.3, 1), 2), source = 1, bcs = dirichlet_bc(0)))
  b <- solve_pde(elliptic_problem(m, diffusion = tensor_coefficient(function(x, y) 2 + 0 * x, 0.3, 1),
                                  source = 1, bcs = dirichlet_bc(0)))
  expect_lt(max(abs(a$u$values - b$u$values)), 1e-13)
  u <- function(x, y) sin(pi * x) * sin(pi * y)
  f <- function(x, y) -pi * cos(pi * x) * sin(pi * y) + (2 + x) * pi^2 * u(x, y) - pi^2 * cos(pi * x) * cos(pi * y)
  K <- tensor_coefficient(function(x, y) 1 + x, 0.5, 1)
  e <- vapply(c(8, 16), function(n) errornorm(solve_pde(elliptic_problem(pde_mesh_rectangle(nx = n), diffusion = K,
                                                        source = f, bcs = dirichlet_bc(0))), u), numeric(1))
  expect_gt(log2(e[1] / e[2]), 1.9)
  expect_error(solve_pde(elliptic_problem(m, diffusion = tensor_coefficient(1, 2, 1), source = 1,
                                          bcs = dirichlet_bc(0))), "positive definite")
  expect_error(estimate_error(b), "constant tensor")
})
