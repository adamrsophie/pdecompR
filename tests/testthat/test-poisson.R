test_that("affine solutions are reproduced to solver precision", {
  ex <- function(x, y) x + 2 * y
  for (p in 1:2) {
    m <- pde_mesh_rectangle(nx = 5, ny = 4, pattern = "alternating")
    s <- solve_pde(poisson_problem(m, source = 0, boundary = dirichlet_bc(ex), degree = p))
    xy <- s$space$dof_coords
    expect_lt(max(abs(s$u$values - ex(xy[, 1], xy[, 2]))), 1e-12)
    expect_lt(errornorm(s, ex, "H1", gradient = function(x, y) cbind(1 + 0 * x, 2 + 0 * y)), 1e-12)
  }
  L <- pde_mesh_polygon(rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2)), h = 0.5)
  s <- solve_pde(poisson_problem(L, source = 0, boundary = dirichlet_bc(ex)))
  expect_lt(max(abs(s$u$values - ex(L$nodes[, 1], L$nodes[, 2]))), 1e-12)
})

test_that("P1 Poisson converges at rates 2 (L2) and 1 (H1)", {
  prob <- poisson_problem(pde_mesh_rectangle(nx = 4), source = f_sin, boundary = dirichlet_bc(0))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_s3_class(cs, "pde_convergence")
  expect_gt(min(cs$rate_L2[-1]), 1.9)
  expect_gt(min(cs$rate_H1_semi[-1]), 0.95)
  expect_lt(cs$L2[3], 2e-3)
})

test_that("P2 Poisson converges at rates 3 (L2) and 2 (H1)", {
  prob <- poisson_problem(pde_mesh_rectangle(nx = 4), source = f_sin, boundary = dirichlet_bc(0), degree = 2)
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(4, 8, 16))
  expect_gt(min(cs$rate_L2[-1]), 2.85)
  expect_gt(min(cs$rate_H1_semi[-1]), 1.9)
})

test_that("finite-difference gradients in errornorm are labeled and accurate", {
  s <- solve_pde(poisson_problem(pde_mesh_rectangle(nx = 8), source = f_sin, boundary = dirichlet_bc(0)))
  e1 <- errornorm(s, u_sin, "H1_semi", gradient = g_sin)
  e2 <- errornorm(s, u_sin, "H1_semi")
  expect_equal(attr(e2, "gradient"), "finite difference")
  expect_equal(as.numeric(e1), as.numeric(e2), tolerance = 1e-7)
  en <- errornorm(s, u_sin, "energy", gradient = g_sin, diffusion = 4)
  expect_equal(as.numeric(en), 2 * as.numeric(e1), tolerance = 1e-12)
})

test_that("mixed Dirichlet and Neumann conditions are enforced", {
  # u = sin(pi x) sin(pi y) with flux k du/dn on the top edge (y = 1):
  # du/dy = pi sin(pi x) cos(pi y) = -pi sin(pi x)
  flux_top <- function(x, y) -pi * sin(pi * x)
  prob <- poisson_problem(pde_mesh_rectangle(nx = 4), source = f_sin,
                          boundary = list(dirichlet_bc(0, c("bottom", "left", "right")),
                                          neumann_bc(flux_top, "top")))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.85)
  expect_gt(min(cs$rate_H1_semi[-1]), 0.95)
})

test_that("variable diffusion converges", {
  k <- function(x, y) 1 + x^2 + y
  # -div(k grad u) with u = sin sin
  f <- function(x, y) {
    ux <- pi * cos(pi * x) * sin(pi * y); uy <- pi * sin(pi * x) * cos(pi * y)
    k(x, y) * 2 * pi^2 * sin(pi * x) * sin(pi * y) - (2 * x * ux + 1 * uy)
  }
  prob <- poisson_problem(pde_mesh_rectangle(nx = 4), diffusion = k, source = f, boundary = dirichlet_bc(0))
  cs <- convergence_study(prob, u_sin, g_sin, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.9)
  expect_error(solve_pde(poisson_problem(prob$mesh, diffusion = function(x, y) x - 0.5, source = 1,
                                         boundary = dirichlet_bc(0))), "positive")
})

test_that("Robin conditions converge", {
  # u = exp(x + y): on x = 1, du/dn = u, so k du/dn + 2 u = 3 u
  ex <- function(x, y) exp(x + y)
  gr <- function(x, y) cbind(exp(x + y), exp(x + y))
  prob <- elliptic_problem(pde_mesh_rectangle(nx = 4), source = function(x, y) -2 * exp(x + y),
                           bcs = list(dirichlet_bc(ex, c("bottom", "top", "left")),
                                      robin_bc(2, function(x, y) 3 * exp(x + y), "right")))
  cs <- convergence_study(prob, ex, gr, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.9)
})

test_that("reaction terms and full Robin problems give SPD systems", {
  ex <- function(x, y) cos(x) * cos(y)
  # -Laplace u + u = 3 u; on the unit square boundary use robin with alpha = 1
  prob <- elliptic_problem(pde_mesh_rectangle(nx = 16), reaction = 1,
                           source = function(x, y) 3 * cos(x) * cos(y),
                           bcs = list(robin_bc(1, function(x, y) {
                             dn <- ifelse(abs(x - 1) < 1e-12, -sin(x) * cos(y),
                                   ifelse(abs(y - 1) < 1e-12, -cos(x) * sin(y),
                                   ifelse(abs(x) < 1e-12, sin(x) * cos(y), cos(x) * sin(y))))
                             dn + cos(x) * cos(y)
                           })))
  s <- solve_pde(prob)
  expect_equal(s$diagnostics$solver$method, "cholesky")
  expect_lt(as.numeric(errornorm(s, ex, "L2")), 1e-3)
})

test_that("pure Neumann problems are detected, checked and gauged", {
  ex <- function(x, y) cos(pi * x) * cos(pi * y)
  m <- pde_mesh_rectangle(nx = 16)
  s <- solve_pde(poisson_problem(m, source = function(x, y) 2 * pi^2 * cos(pi * x) * cos(pi * y)))
  expect_false(is.null(s$diagnostics$gauge))
  expect_equal(s$diagnostics$gauge$type, "zero mean")
  expect_lt(abs(sum(assemble_vector(s$space, 1) * s$u$values)), 1e-12)
  expect_lt(as.numeric(errornorm(s, ex, "L2")), 1e-2)
  expect_error(solve_pde(poisson_problem(m, source = 1)), "incompatible")
  # compatible Neumann data on the boundary
  s2 <- solve_pde(poisson_problem(m, source = 1, boundary = neumann_bc(-0.25)))
  expect_equal(s2$diagnostics$gauge$type, "zero mean")
})

test_that("summary, print and diagnostics report solver information", {
  s <- solve_pde(poisson_problem(pde_mesh_rectangle(nx = 4), source = 1, boundary = dirichlet_bc(0)))
  expect_output(print(s), "pde_solution")
  expect_output(print(summary(s)), "degrees of freedom")
  expect_lt(s$diagnostics$solver$relative_residual, 1e-12)
  expect_equal(s$diagnostics$ndof, 25)
  expect_equal(s$diagnostics$free, 9)
})

test_that("invalid problems are rejected", {
  m <- pde_mesh_rectangle(nx = 2)
  expect_error(poisson_problem(m, boundary = dirichlet_bc(0, "inlet")), "Unknown")
  expect_error(poisson_problem(m, source = "a"), "finite number")
  expect_error(dirichlet_bc(NA), "finite")
  expect_error(robin_bc(-1), "non-negative")
  expect_error(poisson_problem(m, boundary = list(1)), "pde_bc")
})
