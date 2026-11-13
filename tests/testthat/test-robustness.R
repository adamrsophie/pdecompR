test_that("products of two directional derivatives are rejected", {
  V <- fe_space(pde_mesh_rectangle(nx = 3))
  u <- trial_function(V); v <- test_function(V)
  expect_error(integrate_domain(dot(c(1, 0), grad(u)) * dot(c(1, 0), grad(v))), "directional derivatives")
})

test_that("the estimator is exact for a discontinuous coefficient with a P1 solution", {
  k <- function(x, y) ifelse(x < 0.5, 1, 100)
  s <- solve_pde(elliptic_problem(pde_mesh_rectangle(nx = 8), diffusion = k, source = 0,
                                  bcs = list(dirichlet_bc(0, "left"), dirichlet_bc(1, "right"))))
  expect_lt(estimate_error(s)$total, 1e-8)
})

test_that("the estimator lets Dirichlet conditions win on shared edges", {
  m <- pde_mesh_rectangle(nx = 8)
  a <- estimate_error(solve_pde(poisson_problem(m, source = 1,
                                                boundary = list(dirichlet_bc(0), neumann_bc(1, "right")))))
  b <- estimate_error(solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(0))))
  expect_equal(a$total, b$total)
})

test_that("variational problems refuse new meshes and singular systems are detected", {
  m <- pde_mesh_rectangle(nx = 4)
  V <- fe_space(m); u <- trial_function(V); v <- test_function(V)
  vp <- variational_problem(integrate_domain(dot(grad(u), grad(v))), integrate_domain(1 * v))
  expect_error(update_mesh(vp, m), "tied to their function space")
  expect_error(convergence_study(vp, function(x, y) 0, refinements = 1), "variational")
  expect_error(solve_pde(vp), "incompatible")
  expect_error(solve_pde(elliptic_problem(m, source = 1, bcs = robin_bc(0, 0))), "incompatible")
  expect_error(solve_pde(elasticity_problem(m, bcs = list(displacement_bc(0, "left", components = "x"),
                                                          traction_bc(c(1, 0), "right")))),
               "rigid motions")
})

test_that("pure Neumann advection problems use the left null vector", {
  m <- pde_mesh_rectangle(nx = 16)
  # int f exp(-x) = 0 is the solvability condition for -u'' + u' = f with zero flux
  f <- function(x, y) exp(x) * (x - 0.5)
  s <- expect_no_warning(solve_pde(elliptic_problem(m, advection = c(1, 0), source = f)))
  expect_equal(s$diagnostics$gauge$type, "zero mean")
  expect_warning(solve_pde(elliptic_problem(m, advection = c(1, 0), source = function(x, y) x - 0.5)),
                 "compatibility")
})

test_that("time-dependent natural boundary data work in parabolic problems", {
  m <- pde_mesh_rectangle(nx = 4)
  h <- heat_problem(m, initial = 0, time = c(0, 1),
                    bcs = list(dirichlet_bc(0, "left"), robin_bc(1, function(x, y, t) t, "right")))
  expect_gt(max(solve_pde(h, dt = 0.25)$values), 0)
  r <- reaction_diffusion_problem(m, reaction = function(u, x, y, t) 0 * u, initial = list(u = 0),
                                  bcs = list(u = neumann_bc(function(x, y, t) t, "right")), time = c(0, 1))
  expect_gt(max(solve_pde(r, dt = 0.25)$values$u), 0)
})

test_that("extra function arguments are not mistaken for time", {
  m <- pde_mesh_rectangle(nx = 2)
  hp <- heat_problem(m, source = function(x, y, amp = 5) amp + 0 * x, initial = 0, time = c(0, 1))
  expect_equal(mean(solve_pde(hp, dt = 0.5, method = "backward-euler")$values[, 3]), 5, tolerance = 1e-10)
  expect_error(solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(function(x, y, t) t))),
               "stationary")
})

test_that("non-converged linear solves warn", {
  m <- pde_mesh_rectangle(nx = 8)
  expect_warning(solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(0)), method = "cg", maxit = 3),
                 "not converged")
})

test_that("finite-element coefficients on other spaces and reference solutions work", {
  m <- pde_mesh_rectangle(nx = 8)
  kf <- interpolate(function(x, y) 1 + x, fe_space(pde_mesh_rectangle(nx = 16)))
  a <- solve_pde(poisson_problem(m, diffusion = kf, source = 1, boundary = dirichlet_bc(0), degree = 2))
  b <- solve_pde(poisson_problem(m, diffusion = function(x, y) 1 + x, source = 1,
                                 boundary = dirichlet_bc(0), degree = 2))
  expect_lt(max(abs(a$u$values - b$u$values)), 1e-10)
  fine <- solve_pde(poisson_problem(pde_mesh_rectangle(nx = 32), source = 1, boundary = dirichlet_bc(0), degree = 2))
  e1 <- errornorm(solve_pde(poisson_problem(pde_mesh_rectangle(nx = 4), source = 1, boundary = dirichlet_bc(0))), fine$u)
  e2 <- errornorm(solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(0))), fine$u)
  expect_gt(log2(e1 / e2), 1.7)
  kp <- interpolate(function(x, y) 1 + x, fe_space(m))
  expect_true(is.finite(estimate_error(solve_pde(poisson_problem(m, diffusion = kp, source = 1,
                                                                 boundary = dirichlet_bc(0))))$total))
})

test_that("relabeled meshes and marker typos are caught", {
  m <- set_boundary_marker(pde_mesh_rectangle(nx = 4), "inlet", function(x, y) x < 1e-8)
  p <- poisson_problem(m, source = 1, boundary = dirichlet_bc(0, "inlet"))
  expect_error(convergence_study(p, function(x, y) 0, mesh_sizes = c(4, 8)), "relabeled")
  expect_error(boundary_nodes(m, c("boundary", "typo")), "typo")
})

test_that("reaction-diffusion inputs are matched by species", {
  m <- pde_mesh_rectangle(nx = 4)
  rf <- function(u, x, y, t) 0 * u
  expect_error(reaction_diffusion_problem(m, species = c("u", "v"), reaction = rf, initial = list(a = 1, b = 2)),
               "must match")
  expect_error(reaction_diffusion_problem(m, reaction = rf, initial = 0, bcs = list(U = dirichlet_bc(0))),
               "species names")
  r <- reaction_diffusion_problem(m, species = c("a", "b"), reaction = rf, initial = list(a = 1, b = 2),
                                  bcs = list(a = dirichlet_bc(0)))
  expect_output(print(r), "a: dirichlet")
  expect_s3_class(update_mesh(r, pde_mesh_rectangle(nx = 2)), "reaction_diffusion_problem")
  s <- solve_pde(r, dt = 0.5)
  expect_equal(c(s$diagnostics$ndof, s$diagnostics$free), c(50, 34))
  r1 <- reaction_diffusion_problem(m, diffusion = function(x, y) 1 + x, reaction = rf,
                                   initial = interpolate(function(x, y) x, fe_space(m)))
  expect_s3_class(solve_pde(r1, dt = 0.5), "pde_transient")
})

test_that("SUPG accepts tensor diffusion and Crank-Nicolson starts smoothly", {
  m <- pde_mesh_rectangle(nx = 8)
  s <- solve_pde(elliptic_problem(m, diffusion = diag(c(0.01, 0.02)), advection = c(1, 0),
                                  stabilization = "supg", source = 1, bcs = dirichlet_bc(0)))
  expect_true(all(is.finite(s$u$values)))
  h <- heat_problem(pde_mesh_rectangle(nx = 16), initial = 1, bcs = dirichlet_bc(0), time = c(0, 0.5))
  expect_gt(min(solve_pde(h, dt = 0.05)$values), -1e-2)
  expect_lt(min(solve_pde(h, dt = 0.05, startup = 0)$values), -0.5)
})

test_that("summaries, plots and VTK output work for all solution types", {
  m <- pde_mesh_rectangle(nx = 4)
  s <- solve_pde(nonlinear_problem(m, diffusion = function(u) u^2, source = 1, bcs = dirichlet_bc(1)))
  expect_output(print(summary(s)), "newton")
  pf <- function(theta, mesh) poisson_problem(mesh, diffusion = theta[1], source = 1, boundary = dirichlet_bc(0))
  expect_output(print(summary(pde_sensitivity(pf, 1, mesh = m)$solution)), "cholesky")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot(interpolate(1, fe_space(m))))
  expect_silent(plot(interpolate(function(x, y) x * y, fe_space(m, 2))))
  f <- tempfile(fileext = ".vtk")
  write_vtk(solve_pde(elasticity_problem(m, bcs = list(displacement_bc(0, "left"), traction_bc(c(0, -1), "right")))), f)
  expect_true(any(grepl("^VECTORS u double", readLines(f))))
  h <- solve_pde(heat_problem(m, initial = 1, bcs = dirichlet_bc(0), time = c(0, 0.2)), dt = 0.1)
  write_vtk(h, f)
  expect_true(any(grepl("^SCALARS u double", readLines(f))))
  expect_equal(evaluate(s$u, c(0.5, 0.5)), evaluate(s$u, 0.5, 0.5))
})

test_that("weak forms reject grad of coefficient products and mixed spaces", {
  m <- pde_mesh_rectangle(nx = 3)
  V <- fe_space(m); W <- fe_space(m, 2)
  u <- trial_function(V); v <- test_function(V)
  expect_error(grad(function(x, y) 1 + x) , "trial or test")
  expect_error(grad((function(x, y) 1 + x) * u), "spatial coefficient")
  expect_error(dot(grad(u), grad(test_function(W))), "same function space")
  expect_error(variational_problem(integrate_domain(u * v), integrate_domain(1 * test_function(W))),
               "same function space")
})

test_that("invalid data are rejected", {
  m <- pde_mesh_rectangle(nx = 2)
  expect_error(solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(function(x, y) NA_real_ + x))),
               "non-finite")
  expect_error(fe_function(fe_space(m), factor(rep("10", 9))), "numeric")
  expect_error(pde_observations(0.5, 0.5, "a"), "finite numbers")
  expect_error(p_laplace_problem(m, p = 3, epsilon = 0), "epsilon")
  expect_error(pde_mesh(m$nodes, m$triangles, cbind(c(0, 2), c(6, 3))), "node indices")
  sq <- pde_domain("polygon", vertices = cbind(c(0, 1, 1, 0), c(0, 0, 1, 1)) * 1e-8)
  expect_s3_class(sq, "pde_domain")
})

test_that("iterative solvers return zero for a zero right-hand side", {
  A <- Matrix::Diagonal(2)
  for (meth in c("cg", "bicgstab", "gmres")) {
    expect_equal(solve_linear_system(A, c(0, 0), meth, x0 = c(1, 1))$x, c(0, 0))
  }
  expect_error(solve_linear_system(A, c(1, 1), "gmres", restart = 0), "restart")
})

test_that("initial functions are transferred by geometry, constant callbacks recycle", {
  a <- pde_mesh_rectangle(c(0, 2), c(0, 1), nx = 2, ny = 2)
  b <- pde_mesh_rectangle(nx = 2, ny = 2)
  u0 <- interpolate(function(x, y) x, fe_space(a))
  h <- heat_problem(b, initial = u0, time = c(0, 0.1))
  expect_equal(solve_pde(h, dt = 0.1)$values[, 1], b$nodes[, 1])
  s <- solve_pde(nonlinear_problem(pde_mesh_rectangle(nx = 4), diffusion = function(u) 1, source = 1,
                                   bcs = dirichlet_bc(0)))
  expect_true(s$converged)
})

test_that("the estimator sees Dirichlet data the space cannot represent", {
  s <- solve_pde(poisson_problem(pde_mesh_rectangle(nx = 1), source = 0,
                                 boundary = dirichlet_bc(function(x, y) x^2 - y^2)))
  expect_gt(estimate_error(s)$total, 0.1)
})
