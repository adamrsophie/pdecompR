make_system <- function(n = 10, adv = NULL) {
  V <- fe_space(pde_mesh_rectangle(nx = n))
  A <- assemble_matrix(V, "stiffness") + assemble_matrix(V, "mass")
  if (!is.null(adv)) A <- A + assemble_matrix(V, "advection", adv)
  list(A = A, b = assemble_vector(V, function(x, y) 1 + x))
}

test_that("direct and iterative solvers agree on SPD systems", {
  s <- make_system()
  ref <- solve_linear_system(s$A, s$b)
  expect_equal(ref$info$method, "cholesky")
  expect_lt(ref$info$relative_residual, 1e-12)
  for (m in c("lu", "cg", "bicgstab", "gmres")) {
    r <- solve_linear_system(s$A, s$b, method = m, tol = 1e-12)
    expect_equal(r$x, ref$x, tolerance = 1e-9, label = m)
    if (m != "lu") expect_true(r$info$iterations > 0)
  }
})

test_that("nonsymmetric systems use LU or Krylov methods", {
  s <- make_system(adv = c(20, 5))
  r <- solve_linear_system(s$A, s$b)
  expect_equal(r$info$method, "lu")
  expect_false(r$info$symmetric)
  expect_error(solve_linear_system(s$A, s$b, method = "cg"), "symmetric")
  expect_error(solve_linear_system(s$A, s$b, method = "cholesky"), "symmetric")
  g <- solve_linear_system(s$A, s$b, method = "gmres", tol = 1e-12, restart = 30)
  expect_equal(g$x, r$x, tolerance = 1e-8)
  bi <- solve_linear_system(s$A, s$b, method = "bicgstab", tol = 1e-12)
  expect_equal(bi$x, r$x, tolerance = 1e-8)
})

test_that("indefinite symmetric systems fall back from Cholesky", {
  A <- Matrix::Matrix(c(1, 2, 2, 1), 2, sparse = TRUE)
  r <- solve_linear_system(A, c(1, 1))
  expect_equal(r$info$method, "lu")
  expect_match(r$info$note, "not positive definite")
  expect_equal(r$x, c(1, 1) / 3)
  expect_error(solve_linear_system(A, c(1, 1), method = "cholesky"), "positive definite")
})

test_that("iterative solvers report non-convergence", {
  s <- make_system()
  r <- solve_linear_system(s$A, s$b, method = "cg", tol = 1e-14, maxit = 2)
  expect_equal(r$info$status, "not converged")
})

test_that("Poisson solves are solver independent", {
  prob <- poisson_problem(pde_mesh_rectangle(nx = 8), source = f_sin, boundary = dirichlet_bc(0))
  a <- solve_pde(prob)
  b <- solve_pde(prob, method = "cg", tol = 1e-13)
  expect_equal(a$u$values, b$u$values, tolerance = 1e-10)
})

test_that("smoothed aggregation AMG gives mesh-independent CG iterations", {
  L <- pde_domain("polygon", vertices = cbind(c(0, 2, 2, 1, 1, 0), c(0, 0, 1, 1, 2, 2)))
  it <- vapply(c(0.1, 0.05), function(h) {
    p <- poisson_problem(pde_mesh_polygon(L, h = h), source = 1, boundary = dirichlet_bc(0), degree = 2)
    s <- solve_pde(p, method = "cg", preconditioner = amg_preconditioner(), tol = 1e-10)
    d <- solve_pde(p)
    expect_equal(s$u$values, d$u$values, tolerance = 1e-7)
    s$diagnostics$solver$iterations
  }, numeric(1))
  expect_lt(it[2], 1.5 * it[1])
  expect_lt(it[2], 60)
  A <- assemble_matrix(fe_space(pde_mesh_rectangle(nx = 4)), "stiffness") + assemble_matrix(fe_space(pde_mesh_rectangle(nx = 4)), "mass")
  s <- solve_linear_system(A, rep(1, nrow(A)), method = "cg", preconditioner = amg_preconditioner())
  expect_equal(s$info$status, "converged")
  expect_error(amg_preconditioner(strength = 1), "strength")
})
