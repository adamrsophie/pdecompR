test_that("Taylor-Hood Stokes converges at the optimal rates", {
  a <- function(x) x^2 - 2 * x^3 + x^4; a1 <- function(x) 2 * x - 6 * x^2 + 4 * x^3
  a2 <- function(x) 2 - 12 * x + 12 * x^2; a3 <- function(x) -12 + 24 * x
  nu <- 0.5
  ux <- function(x, y) a(x) * a1(y)
  pex <- function(x, y) x * y - 0.25
  f <- function(x, y) cbind(-nu * (a2(x) * a1(y) + a(x) * a3(y)) + y, nu * (a3(x) * a(y) + a1(x) * a2(y)) + x)
  r <- vapply(c(4, 8), function(n) {
    s <- solve_pde(stokes_problem(pde_mesh_rectangle(nx = n), viscosity = nu, body_force = f, bcs = velocity_bc(c(0, 0))))
    expect_equal(s$diagnostics$pressure_gauge, "zero mean")
    c(errornorm(s$velocity$ux, ux), errornorm(s$pressure, pex))
  }, numeric(2))
  expect_gt(log2(r[1, 1] / r[1, 2]), 2.8)
  expect_gt(log2(r[2, 1] / r[2, 2]), 1.9)
  its <- vapply(c(4, 8, 16), function(n) {
    s <- solve_pde(stokes_problem(pde_mesh_rectangle(nx = n), viscosity = nu, body_force = f, bcs = velocity_bc(c(0, 0))),
                   method = "gmres")
    expect_equal(s$diagnostics$solver$status, "converged")
    s$diagnostics$solver$iterations
  }, numeric(1))
  expect_lt(max(its) - min(its), 15)
})

test_that("Poiseuille flow with do-nothing outflow is reproduced exactly", {
  m <- pde_mesh_rectangle(c(0, 2), c(0, 1), nx = 4, ny = 2)
  s <- solve_pde(stokes_problem(m, bcs = list(velocity_bc(c(0, 0), c("top", "bottom")),
                                              velocity_bc(function(x, y) cbind(y * (1 - y), 0), "left"))))
  expect_lt(errornorm(s$velocity$ux, function(x, y) y * (1 - y)), 1e-12)
  expect_lt(errornorm(s$pressure, function(x, y) -2 * x + 4), 1e-11)
  expect_lt(s$diagnostics$divergence_norm, 1e-12)
  expect_output(print(s), "Taylor-Hood")
  expect_error(stokes_problem(m, viscosity = -1), "positive")
})

test_that("closed Stokes problems need zero net boundary flux", {
  m <- pde_mesh_rectangle(nx = 4)
  expect_error(solve_pde(stokes_problem(m, bcs = list(velocity_bc(c(0, 0)), velocity_bc(c(1, 1), "top")))),
               "net boundary flux")
  expect_s3_class(solve_pde(stokes_problem(m, bcs = list(velocity_bc(c(0, 0)), velocity_bc(c(1, 0), "top")))),
                  "pde_stokes_solution")
})
