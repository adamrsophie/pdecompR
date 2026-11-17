# Manufactured displacement with lambda = 1.5, mu = 1 on the unit square.
el_lam <- 1.5; el_mu <- 1
el_E <- el_mu * (3 * el_lam + 2 * el_mu) / (el_lam + el_mu); el_nu <- el_lam / (2 * (el_lam + el_mu))
el_exact <- function(x, y) cbind(sin(pi * x) * sin(pi * y), x * y * (1 - x) * (1 - y))
el_grad <- list(function(x, y) cbind(pi * cos(pi * x) * sin(pi * y), pi * sin(pi * x) * cos(pi * y)),
                function(x, y) cbind(y * (1 - y) * (1 - 2 * x), x * (1 - x) * (1 - 2 * y)))
el_force <- function(x, y) {
  uxx <- -pi^2 * sin(pi * x) * sin(pi * y); uxy <- pi^2 * cos(pi * x) * cos(pi * y)
  vxx <- -2 * y * (1 - y); vyy <- -2 * x * (1 - x); vxy <- (1 - 2 * x) * (1 - 2 * y)
  cbind(-(el_mu * 2 * uxx + (el_lam + el_mu) * (uxx + vxy)),
        -(el_mu * (vxx + vyy) + (el_lam + el_mu) * (uxy + vyy)))
}

test_that("vector P1 and P2 elasticity converge at optimal rates", {
  p <- elasticity_problem(pde_mesh_rectangle(nx = 4), E = el_E, nu = el_nu, body_force = el_force,
                          bcs = displacement_bc(c(0, 0)))
  cs <- convergence_study(p, el_exact, el_grad, mesh_sizes = c(8, 16, 32))
  expect_gt(min(cs$rate_L2[-1]), 1.9)
  expect_gt(min(cs$rate_H1_semi[-1]), 0.95)
  expect_equal(cs$ndof[1], 2 * 81)
  p2 <- elasticity_problem(pde_mesh_rectangle(nx = 4), E = el_E, nu = el_nu, body_force = el_force,
                           bcs = displacement_bc(c(0, 0)), degree = 2)
  cs2 <- convergence_study(p2, el_exact, el_grad, mesh_sizes = c(8, 16))
  expect_gt(cs2$rate_L2[2], 2.85)
  expect_gt(cs2$rate_H1_semi[2], 1.9)
})

test_that("affine displacements are reproduced exactly (patch test)", {
  m <- pde_mesh_polygon(rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2)), h = 0.4)
  lin <- function(x, y) cbind(0.1 * x - 0.3 * y + 1, 0.2 * x + 0.05 * y)
  for (deg in 1:2) {
    s <- solve_pde(elasticity_problem(m, E = 3, nu = 0.25, bcs = displacement_bc(lin), degree = deg))
    xy <- s$space$dof_coords
    expect_lt(max(abs(cbind(s$ux$values, s$uy$values) - lin(xy[, 1], xy[, 2]))), 1e-12)
  }
})

test_that("a P2 cantilever matches Timoshenko beam theory", {
  m <- pde_mesh_rectangle(c(0, 10), c(-0.5, 0.5), nx = 20, ny = 4)
  s <- solve_pde(elasticity_problem(m, E = 1000, nu = 0.3, plane = "stress", degree = 2,
                                    bcs = list(displacement_bc(c(0, 0), "left"),
                                               traction_bc(c(0, -1), "right"))))
  # bending P L^3 / (3 E I) plus shear P L / (kappa G A)
  theory <- 1 * 10^3 / (3 * 1000 / 12) + 10 / (5 / 6 * 1000 / 2.6)
  expect_equal(-evaluate(s$uy, 10, 0), theory, tolerance = 0.01)
  expect_equal(s$diagnostics$solver$method, "cholesky")
  expect_output(print(s), "displacement")
  grDevices::pdf(NULL); withr::defer(grDevices::dev.off())
  expect_silent(plot(s))
})

test_that("rollers, rigid motions and parameters are checked", {
  m <- pde_mesh_rectangle(nx = 4)
  # uniaxial tension with rollers: exact u_x = x * sigma / E, plane stress
  s <- solve_pde(elasticity_problem(m, E = 2, nu = 0.25, plane = "stress",
                                    bcs = list(displacement_bc(0, "left", components = "x"),
                                               displacement_bc(0, "bottom", components = "y"),
                                               traction_bc(c(1, 0), "right"))))
  xy <- s$space$dof_coords
  expect_lt(max(abs(s$ux$values - xy[, 1] / 2)), 1e-12)
  expect_lt(max(abs(s$uy$values + 0.25 * xy[, 2] / 2)), 1e-12)
  expect_error(solve_pde(elasticity_problem(m, bcs = traction_bc(c(1, 0), "right"))), "rigid motions")
  expect_error(solve_pde(elasticity_problem(m, nu = 0.5, bcs = displacement_bc(c(0, 0)))), "Poisson")
  expect_error(solve_pde(elasticity_problem(m, E = -1, bcs = displacement_bc(c(0, 0)))), "positive")
  expect_error(displacement_bc(c(1, 2, 3)), "length-2")
  expect_error(elasticity_problem(m, bcs = dirichlet_bc(0)), "displacement_bc")
  expect_error(elasticity_problem(m, bcs = displacement_bc(0, "nowhere")), "nowhere")
  expect_output(print(displacement_bc(0, "left", "x")), "displacement \\[x\\]")
})

test_that("the mixed formulation does not lock as nu approaches 1/2", {
  ps <- quote((x * (1 - x) * y * (1 - y))^2)
  ux_e <- D(ps, "y"); uy_e <- call("-", D(ps, "x"))
  lap <- function(e) call("+", D(D(e, "x"), "x"), D(D(e, "y"), "y"))
  lx <- lap(ux_e); ly <- lap(uy_e)
  err <- function(nu, fm, deg) {
    mu <- 1 / (2 * (1 + nu))
    s <- solve_pde(elasticity_problem(pde_mesh_rectangle(nx = 8), nu = nu, degree = deg, formulation = fm,
                                      body_force = function(x, y) cbind(-mu * eval(lx), -mu * eval(ly)),
                                      bcs = displacement_bc(0)))
    xy <- s$space$dof_coords
    max(abs(s$ux$values - eval(ux_e, list(x = xy[, 1], y = xy[, 2]))))
  }
  m2 <- err(0.49, "mixed", 2); m4 <- err(0.4999, "mixed", 2)
  expect_lt(abs(m4 / m2 - 1), 0.05)
  expect_lt(m4, 1.5 * err(0.3, "mixed", 2))
  expect_gt(err(0.4999, "displacement", 1), 50 * m4)
  expect_error(elasticity_problem(pde_mesh_rectangle(nx = 2), formulation = "mixed"), "degree = 2")
  ep <- elasticity_problem(pde_mesh_rectangle(nx = 4), nu = 0.49, degree = 2, formulation = "mixed",
                           body_force = c(0, -1), bcs = displacement_bc(0, "left"))
  a <- solve_pde(ep); g <- solve_pde(ep, method = "gmres", tol = 1e-12)
  expect_equal(g$ux$values, a$ux$values, tolerance = 1e-8)
  expect_s3_class(a$pressure, "fe_function")
  expect_true(all(is.finite(elasticity_stress(a)$von_mises)))
  expect_error(solve_eigen(ep), "displacement formulation")
  expect_error(elastodynamics_problem(ep), "displacement formulation")
})
