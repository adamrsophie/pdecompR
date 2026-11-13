test_that("the weak-form language reproduces the Poisson solver", {
  m <- pde_mesh_rectangle(nx = 8)
  V <- fe_space(m)
  u <- trial_function(V); v <- test_function(V)
  a <- integrate_domain(dot(grad(u), grad(v)))
  L <- integrate_domain(f_sin * v)
  s1 <- solve_pde(variational_problem(a, L, dirichlet_bc(0)))
  s2 <- solve_pde(poisson_problem(m, source = f_sin, boundary = dirichlet_bc(0)))
  expect_equal(s1$u$values, s2$u$values, tolerance = 1e-12)
})

test_that("forms with reaction, advection and boundary terms match the elliptic solver", {
  m <- pde_mesh_rectangle(nx = 6)
  V <- fe_space(m)
  u <- trial_function(V); v <- test_function(V)
  k <- function(x, y) 1 + x
  b <- c(1, -0.5)
  a <- integrate_domain(k * dot(grad(u), grad(v)) + 2 * u * v + dot(b, grad(u)) * v) +
    integrate_boundary(3 * u * v, "right")
  L <- integrate_domain(1 * v) + integrate_boundary((function(x, y) y) * v, "right") +
    integrate_boundary(0.5 * v, "top")
  s1 <- solve_pde(variational_problem(a, L, dirichlet_bc(0, "left")))
  prob <- elliptic_problem(m, diffusion = k, reaction = 2, advection = b, source = 1,
                           bcs = list(dirichlet_bc(0, "left"), robin_bc(3, function(x, y) y, "right"),
                                      neumann_bc(0.5, "top")))
  s2 <- solve_pde(prob)
  expect_equal(s1$u$values, s2$u$values, tolerance = 1e-12)
})

test_that("algebra on expressions and forms works", {
  V <- fe_space(pde_mesh_rectangle(nx = 3))
  u <- trial_function(V); v <- test_function(V)
  A1 <- assemble_form(integrate_domain(dot(grad(u), grad(v)) - u * v), "bilinear")
  A2 <- assemble_matrix(V, "stiffness") - assemble_matrix(V, "mass")
  expect_equal(as.matrix(A1), as.matrix(A2), tolerance = 1e-13)
  A3 <- assemble_form(integrate_domain(dot(grad(u), grad(v))) - integrate_domain(u * v), "bilinear")
  expect_equal(as.matrix(A3), as.matrix(A2), tolerance = 1e-13)
  A4 <- assemble_form(integrate_domain(u * dot(c(1, 2), grad(v))), "bilinear")
  expect_equal(as.matrix(A4), as.matrix(Matrix::t(assemble_matrix(V, "advection", c(1, 2)))), tolerance = 1e-13)
  b <- assemble_form(integrate_domain(dot(vector_coefficient(function(x, y) cbind(x, y)), grad(v))), "linear")
  w <- interpolate(function(x, y) x, V)$values
  expect_equal(sum(b * w), 0.5, tolerance = 1e-13)
})

test_that("ill-formed forms are rejected", {
  V <- fe_space(pde_mesh_rectangle(nx = 2))
  u <- trial_function(V); v <- test_function(V)
  expect_error(u * u, "twice")
  expect_error(grad(u) * grad(v), "dot")
  expect_error(grad(grad(u)), "u or v")
  expect_error(assemble_form(integrate_domain(u * v), "linear"), "bilinear term")
  expect_error(integrate_boundary(u * v, "nope"), "Unknown")
  expect_error(assemble_form(integrate_boundary(dot(grad(u), grad(v))), "bilinear"), "gradients")
  expect_error(variational_problem(integrate_domain(u * v), integrate_domain(v), neumann_bc(1)), "Dirichlet")
  expect_error(grad(1), "trial or test")
})
