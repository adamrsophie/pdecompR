# Numerical contracts: compatibility, nullspaces, Krylov status, nested errors,
# zero-unknown systems, Newton trial points and input validation.


test_that('positive forcing is rejected for incompatible advection Neumann data at any domain scale', {
  for (L in c(1, 100)) {
    m <- pde_mesh_rectangle(x=c(0,L), y=c(0,L), nx=16)
    # Coefficients transform consistently from the unit-square PDE:
    # u(X,Y)=v(X/L,Y/L), k=L^2, b=(L,0), f=1.
    p <- elliptic_problem(m, diffusion=L^2, advection=c(L,0), source=1)
    expect_error(solve_pde(p), 'incompatible|compatibility|no solution')
  }
})

test_that('small positive reaction retains the unique constant solution', {
  p <- elliptic_problem(pde_mesh_rectangle(nx=4),
                        diffusion=1, reaction=1e-9, source=1e-9)
  # The exact solution is u=1. A loose tolerance permits conditioning error;
  # it does not permit treating this uniquely solvable problem as singular.
  s <- solve_pde(p)
  expect_lt(max(abs(s$u$values-1)), 1e-3)
  expect_null(s$diagnostics$gauge)
})

test_that('a nonfinite prescribed displacement is never silently left unconstrained', {
  m <- pde_mesh_rectangle(nx=2)
  for (invalid in c(NA_real_, NaN, Inf, -Inf)) {
    bad <- local({z <- invalid; function(x,y) cbind(ifelse(x==0 & y==0,z,x),y)})
    expect_error(solve_pde(elasticity_problem(m,bcs=displacement_bc(bad))))
  }
})

test_that('fully prescribed elasticity returns the imposed displacement', {
  p <- elasticity_problem(pde_mesh_rectangle(nx=1),
                          bcs=displacement_bc(function(x,y) cbind(x,y)))
  s <- solve_pde(p)
  xy <- s$space$dof_coords
  expect_equal(s$ux$values,xy[,1],tolerance=1e-12,ignore_attr=TRUE)
  expect_equal(s$uy$values,xy[,2],tolerance=1e-12,ignore_attr=TRUE)
  expect_equal(s$diagnostics$free,0)
})

test_that('Krylov methods preserve identity solutions at extreme representable scales', {
  for (meth in c('cg','bicgstab','gmres')) {
    for (mag in c(1e-200,1e200)) {
      b <- mag*c(1,2)
      out <- tryCatch(solve_linear_system(Matrix::Diagonal(2),b,method=meth),
                      error=function(e)e)
      # A documented numerical-range rejection is preferable to false success.
      # Supporting this well-conditioned identity system is the desired fix.
      if (inherits(out,'error')) {
        expect_true(nzchar(conditionMessage(out)))
      } else {
        expect_true(all(is.finite(out$x)))
        expect_lt(max(abs(out$x/mag-c(1,2))),1e-9)
        if (out$info$status %in% c('converged','solved')) {
          expect_true(is.finite(out$info$relative_residual))
          expect_lte(out$info$relative_residual,1e-10)
        }
      }
    }
  }
})

test_that('Krylov convergence status agrees with the explicit residual', {
  n <- 8
  A <- Matrix::Matrix(1/outer(seq_len(n),seq_len(n),function(i,j) i+j-1),sparse=TRUE)
  b <- sin(seq_len(n)); tol <- 1e-10
  for (meth in c('cg','bicgstab')) {
    out <- solve_linear_system(A,b,method=meth,tol=tol,maxit=2000)
    true_rel <- sqrt(sum((b-as.numeric(A %*% out$x))^2))/sqrt(sum(b^2))
    expect_equal(out$info$relative_residual,true_rel,tolerance=1e-12)
    if (out$info$status %in% c('converged','solved')) {
      expect_true(is.finite(true_rel))
      expect_lte(true_rel,tol*(1+1e-6))
    } else {
      # Honest non-convergence on this ill-conditioned problem is acceptable.
      expect_true(nzchar(out$info$status))
    }
  }
})


test_that("fully prescribed heat returns the time-dependent data", {
 m <- pde_mesh_rectangle(nx=1)
 p <- heat_problem(m,initial=0,bcs=dirichlet_bc(function(x,y,t)t),time=c(0,1))
 s <- solve_pde(p,dt=.5)
 expect_equal(as.numeric(s$values[,ncol(s$values)]),rep(1,4),tolerance=1e-12)
})

test_that("fully prescribed IMEX reaction diffusion has no linear solve", {
 m <- pde_mesh_rectangle(nx=1)
 p <- reaction_diffusion_problem(m,reaction=function(u,x,y,t)0*u,
   initial=list(u=0),bcs=list(u=dirichlet_bc(0)),time=c(0,1))
 s <- solve_pde(p,dt=.5)
 expect_equal(as.numeric(s$values$u),rep(0,length(s$values$u)),tolerance=1e-12)
})

test_that("fully prescribed forward sensitivity differentiates boundary data", {
 m <- pde_mesh_rectangle(nx=1)
 pf <- function(theta) poisson_problem(m,boundary=dirichlet_bc(theta[1]))
 s <- pde_sensitivity(pf,theta=1)
 expect_equal(as.numeric(s$sensitivities[,1]),rep(1,4),tolerance=1e-8)
})

test_that("fully prescribed optimal control minimizes its control cost", {
 m <- pde_mesh_rectangle(nx=1)
 s <- optimal_control(m,target=1,alpha=1)
 expect_equal(s$control$values,rep(0,4),tolerance=1e-12)
 expect_equal(s$state$values,rep(0,4),tolerance=1e-12)
 expect_equal(s$objective,.5,tolerance=1e-12)
})


test_that('nested directional contractions are rejected or preserve direction', {
  V <- fe_space(pde_mesh_rectangle(nx = 2))
  u <- trial_function(V); v <- test_function(V)
  # Rejection is an acceptable current-release fix because this form is not
  # in the documented supported subset. A future implementation may pass
  # the exact directional-energy test instead.
  ans <- tryCatch({
    a <- integrate_domain(dot(c(1, 0) * dot(c(1, 0), grad(u)), grad(v)))
    discretize(variational_problem(a, integrate_domain(0 * v)))$A
  }, error = identity)
  if (inherits(ans, 'error')) {
    expect_match(conditionMessage(ans), 'supported|directional|contraction',
                 ignore.case = TRUE)
  } else {
    nodal_y <- V$dof_coords[, 2]
    expect_equal(as.numeric(crossprod(nodal_y, ans %*% nodal_y)), 0,
                 tolerance = 1e-12)
  }
})

test_that('FE-reference errors do not miss fine-grid support', {
  Vc <- fe_space(pde_mesh_rectangle(nx = 1))
  Vf <- fe_space(pde_mesh_rectangle(nx = 32))
  z <- fe_function(Vc, numeric(Vc$ndof))
  values <- numeric(Vf$ndof)
  center <- which(Vf$dof_coords[, 1] == 0.5 & Vf$dof_coords[, 2] == 0.5)
  stopifnot(length(center) == 1L)
  values[center] <- 1
  ref <- fe_function(Vf, values)
  expect_equal(as.numeric(errornorm(z, ref, 'L2')), sqrt(1 / 2048),
               tolerance = 1e-12)
  expect_equal(as.numeric(errornorm(z, ref, 'H1_semi')), 2,
               tolerance = 1e-12)
})

test_that('interface estimator does not change under length-unit rescaling', {
  L <- 1e-6
  mesh <- pde_mesh_rectangle(c(-L, L), c(0, 2 * L), nx = 4)
  exact <- function(x, y) ifelse(x < 0, (x + L) * 100 / (101 * L),
                                 100 / 101 + x / (101 * L))
  diffusion <- function(x, y) ifelse(x < 0, 1, 100)
  solution <- solve_pde(poisson_problem(mesh, diffusion = diffusion,
                                        source = 0,
                                        boundary = dirichlet_bc(exact)))
  expect_lt(max(abs(solution$u$values - exact(mesh$nodes[, 1], mesh$nodes[, 2]))),
            1e-12)
  expect_lt(estimate_error(solution)$total, 1e-7)
})



expect_rejected_or_gauged <- function(expr) {
  result <- tryCatch(force(expr), error = identity)
  if (inherits(result, "error")) {
    testthat::expect_match(conditionMessage(result),
                          "[Nn]eumann|[Nn]ull|[Ss]ingular|[Uu]nique|[Gg]auge")
  } else {
    testthat::expect_false(is.null(result$diagnostics$gauge),
                          info = "A nonunique solution must report an explicit gauge.")
  }
}

testthat::test_that("Newton backtracks after exponential overflow at a trial point", {
  mesh <- pde_mesh_rectangle(nx = 3)
  problem <- nonlinear_problem(mesh, reaction = exp,
    reaction_derivative = exp, source = 1000, initial = 0)
  solution <- solve_pde(problem)
  testthat::expect_true(solution$converged)
  testthat::expect_lt(max(abs(solution$u$values - log(1000))), 1e-8)
})

testthat::test_that("parameter estimation supports a fixed parameter and a free one", {
  mesh <- pde_mesh_rectangle(nx = 3)
  problem_fun <- function(theta) poisson_problem(mesh,
    diffusion = theta[1], source = theta[2], boundary = dirichlet_bc(0))
  truth <- solve_pde(problem_fun(c(1, 2)))
  x <- c(0.25, 0.5, 0.75)
  y <- rep(0.5, length(x))
  observations <- pde_observations(x, y, evaluate(truth$u, x, y), sd = 0.01)
  estimate <- estimate_parameters(problem_fun, c(1, 1), observations,
    lower = c(1, 0), upper = c(1, 10),
    control = list(factr = 1e5, pgtol = 1e-10))
  testthat::expect_equal(estimate$estimate[1], 1, tolerance = 1e-12)
  testthat::expect_equal(estimate$estimate[2], 2, tolerance = 1e-5)
  testthat::expect_equal(estimate$convergence, 0L)
})

testthat::test_that("zero Robin alpha does not bypass nonlinear nullspace handling", {
  mesh <- pde_mesh_rectangle(nx = 3)
  problem <- nonlinear_problem(mesh, diffusion = 1, source = 0,
    bcs = robin_bc(alpha = 0), initial = 7)
  expect_rejected_or_gauged(solve_pde(problem))
})

testthat::test_that("identically zero reaction does not bypass nonlinear nullspace handling", {
  mesh <- pde_mesh_rectangle(nx = 3)
  problem <- nonlinear_problem(mesh, diffusion = 1,
    reaction = function(u) 0 * u, source = 0, initial = 7)
  expect_rejected_or_gauged(solve_pde(problem))
})

testthat::test_that("unsupported time-dependent advection is rejected when manufactured", {
  testthat::skip_if_not_installed("symcompR")
  # Minimal contract-preserving correction: reject time-dependent advection
  # consistently with time-dependent diffusion/reaction in this helper.
  # Replace this test with callback/value tests if time-dependent coefficients
  # are deliberately added as a documented manufactured_solution capability.
  testthat::expect_error(
    manufactured_solution("x", advection = c("t", "0")),
    "[Tt]ime|[Dd]epend|[Aa]dvection")
})

testthat::test_that("reaction diffusion rejects negative evaluated diffusion", {
  mesh <- pde_mesh_rectangle(nx = 3)
  testthat::expect_error({
    problem <- reaction_diffusion_problem(mesh,
      diffusion = function(x, y) rep(-1, length(x)),
      reaction = function(u, x, y, t) 0 * u,
      initial = list(u = 1), time = c(0, 0.001))
    solve_pde(problem, dt = 0.001)
  }, "[Nn]onnegative|[Nn]on-negative|[Dd]iffusion")
})

testthat::test_that("assembly enforces the documented symmetric tensor contract", {
  space <- fe_space(pde_mesh_rectangle(nx = 3))
  diffusion <- rbind(c(1, 4), c(-4, 1))
  # The symmetric part is positive definite, but the supplied tensor is not
  # symmetric. Accepting it would violate the current documented contract.
  testthat::expect_error(assemble_matrix(space, "stiffness", coef = diffusion),
                        "[Ss]ymmetr")
})

test_that("compatible advection Neumann data solve identically at every length scale", {
  out <- sapply(c(1, 100, 1e-4), function(L) {
    m <- pde_mesh_rectangle(c(0, L), c(0, L), nx = 8)
    f <- function(x, y) pi^2 * cos(pi * x / L) - pi * sin(pi * x / L)
    s <- solve_pde(elliptic_problem(m, diffusion = L^2, advection = c(L, 0), source = f))
    c(s$diagnostics$gauge$relative_compatibility, s$diagnostics$solver$relative_residual)
  })
  expect_equal(out[, 2], out[, 1], tolerance = 1e-6)
  expect_equal(out[, 3], out[, 1], tolerance = 1e-6)
})

test_that("disconnected meshes are rejected for pure Neumann and checked per component in elasticity", {
  m <- pde_mesh_rectangle(nx = 2)
  m2 <- pde_mesh(rbind(m$nodes, m$nodes + 2), rbind(m$triangles, m$triangles + nrow(m$nodes)))
  expect_equal(mesh_components(m2)$n, 2L)
  # One zero-mean gauge per component: u = x on each piece, compatible data.
  s <- solve_pde(elliptic_problem(m2, source = 0,
                                  bcs = list(neumann_bc(function(x, y) ifelse(x %in% c(0, 2), -1, ifelse(x %in% c(1, 3), 1, 0))))))
  expect_equal(s$diagnostics$gauge$components, 2L)
  dc <- pdecompR:::dof_components(s$space)
  ex <- s$space$dof_coords[, 1] - c(0.5, 2.5)[dc]
  expect_lt(max(abs(s$u$values - ex)), 1e-10)
  expect_error(solve_pde(elliptic_problem(m2, source = function(x, y) as.numeric(x > 1.5))),
               "component 2")
  m3 <- set_boundary_marker(m2, "first", function(x, y) x <= 1)
  expect_error(solve_pde(elasticity_problem(m3, bcs = displacement_bc(0, "first"))), "component")
})

test_that("finite-element references on unrelated meshes are integrated on the overlay", {
  f <- function(x, y) x^2 * y + sin(3 * x)
  a <- interpolate(f, fe_space(pde_mesh_rectangle(nx = 4)))
  b <- interpolate(f, fe_space(pde_mesh_rectangle(nx = 5, pattern = "left"), 2))
  for (ty in c("L2", "H1_semi")) expect_equal(as.numeric(errornorm(a, b, ty)), as.numeric(errornorm(b, a, ty)),
                                              tolerance = 1e-12)
  sq <- pde_domain("polygon", vertices = cbind(c(0, 1, 1, 0), c(0, 0, 1, 1)))
  g <- function(x, y) 2 * x + y
  expect_lt(errornorm(interpolate(g, fe_space(pde_mesh_rectangle(nx = 3))),
                      interpolate(g, fe_space(pde_mesh_polygon(sq, h = 0.3))), "H1"), 1e-12)
  fine <- interpolate(f, fe_space(refine_uniform(pde_mesh_rectangle(nx = 4))))
  expect_equal(sqrt(sum(pdecompR:::overlay_error(a, fine, "L2", 1, 6))), as.numeric(errornorm(a, fine)),
               tolerance = 1e-10)
  expect_error(errornorm(a, interpolate(f, fe_space(pde_mesh_rectangle(c(0, 2), nx = 3)))), "same domain")
})

test_that("shared validators reject malformed input", {
  m <- pde_mesh_rectangle(nx = 2)
  expect_error(refine_mesh(m, 1.5), "whole-number")
  expect_error(refine_mesh(m, c(TRUE, NA, rep(FALSE, 6))), "TRUE/FALSE")
  expect_error(evaluate(interpolate(1, fe_space(m)), c(0.5, NA), c(0.5, 0.5)), "finite")
  expect_error(evaluation_matrix(fe_space(m), c(0.1, 0.2), 0.3), "same length")
  expect_error(set_boundary_marker(m, "", function(x, y) x == 0), "nonempty")
  expect_error(set_boundary_marker(m, "a", function(x, y) rep(NA, length(x))), "no NA")
  expect_error(solve_pde(heat_problem(m, time = c(0, 1)), dt = 0.5, output_times = numeric(0)), "nonempty")
  expect_error(solve_pde(nonlinear_problem(m, bcs = dirichlet_bc(0)), maxit = 0), "maxit")
  expect_error(mark_elements(c(1, NA)), "finite")
})

test_that("direct solvers are insensitive to diagonal scaling", {
  s <- solve_linear_system(Matrix::Diagonal(x = c(1e-14, 1)), c(1e-14, 1))
  expect_equal(s$info$status, "solved")
  expect_equal(s$x, c(1, 1), tolerance = 1e-14)
})

test_that("variational problems declare the constant nullspace from the form structure", {
  V <- fe_space(pde_mesh_rectangle(nx = 16))
  u <- trial_function(V); v <- test_function(V)
  s <- solve_pde(variational_problem(integrate_domain(dot(grad(u), grad(v)) + 1e-9 * u * v),
                                     integrate_domain(1e-9 * v)))
  expect_null(s$diagnostics$gauge)
  expect_lt(max(abs(s$u$values - 1)), 1e-4)
  expect_error(solve_pde(variational_problem(integrate_domain(dot(grad(u), grad(v))),
                                             integrate_domain(1 * v))), "incompatible")
  s2 <- solve_pde(variational_problem(integrate_domain(dot(grad(u), grad(v)) + 0 * u * v),
                                      integrate_domain(0 * v), nullspace = "constant"))
  expect_false(is.null(s2$diagnostics$gauge))
})

test_that("exact and affine derivatives match finite differences and pass Taylor tests", {
  mesh <- pde_mesh_rectangle(nx = 6)
  pf <- function(th) poisson_problem(mesh, diffusion = th[1], source = th[2], boundary = dirichlet_bc(th[3]))
  ob <- misfit_objective(fe_space(mesh), pde_observations(c(0.3, 0.5), c(0.5, 0.6), c(0.1, 0.2)))
  th <- c(1.5, 2, 0.3)
  dv <- affine_derivative(pf, c(1, 1, 1))
  g_fd <- pde_sensitivity(pf, th, objective = ob, method = "adjoint")$gradient
  g_ad <- pde_sensitivity(pf, th, objective = ob, method = "adjoint", derivative = dv)$gradient
  g_fw <- pde_sensitivity(pf, th, objective = ob, derivative = dv)$gradient
  expect_equal(g_ad, g_fd, tolerance = 1e-7)
  expect_equal(g_fw, g_ad, tolerance = 1e-10)
  tt <- taylor_test(pf, th, ob, derivative = dv)
  expect_gt(min(tt$rate, na.rm = TRUE), 1.9)
  est <- estimate_parameters(function(t) pf(c(t, 2, 0)), 1, pde_observations(c(0.3, 0.5), c(0.5, 0.6),
                             evaluate(solve_pde(pf(c(0.7, 2, 0)))$u, c(0.3, 0.5), c(0.5, 0.6)), sd = 0.01),
                             lower = 0.1)
  expect_equal(est$estimate, 0.7, tolerance = 1e-5)
  expect_true(is.finite(est$identifiability$condition))
})

test_that("reference-element kernels agree with the quadrature loop", {
  L <- pde_domain("polygon", vertices = cbind(c(0, 1, 1, 0.5, 0.5, 0), c(0, 0, 0.5, 0.5, 1, 1)))
  for (deg in 1:2) {
    geo <- pdecompR:::element_geometry(fe_space(pde_mesh_polygon(L, h = 0.2), deg))
    loop <- function(...) { g <- geo; g$gxi <- NULL; pdecompR:::local_bilinear(g, ...) }
    ce <- function(v) matrix(v, geo$m, geo$nq)
    rc <- ce(rep(seq_len(geo$m) / geo$m, geo$nq))
    cases <- list(list(diff = ce(2.5)), list(diff = lapply(list(xx = 2, xy = 0.3, yx = 0.3, yy = 1), ce)),
                  list(mass = rc), list(adv = list(x = ce(1), y = rc)), list(advt = list(x = rc, y = ce(-2))),
                  list(mass = ce(1), diff = geo$x))
    for (a in cases) {
      expect_lt(max(abs(do.call(pdecompR:::local_bilinear, c(list(geo), a)) - do.call(loop, a))), 1e-12)
    }
  }
})

test_that("preconditioners, multigrid and block right-hand sides", {
  ms <- list(pde_mesh_rectangle(nx = 4))
  for (k in 1:3) ms[[k + 1]] <- refine_uniform(ms[[k]])
  its <- vapply(2:4, function(L) {
    prob <- poisson_problem(ms[[L]], source = 1, boundary = dirichlet_bc(0))
    s <- solve_pde(prob, method = "cg", preconditioner = multigrid_preconditioner(ms[1:L]))
    expect_equal(s$diagnostics$solver$status, "converged")
    expect_lt(max(abs(s$u$values - solve_pde(prob)$u$values)), 1e-10)
    s$diagnostics$solver$iterations
  }, numeric(1))
  expect_lte(max(its) - min(its), 3)
  expect_error(multigrid_preconditioner(list(ms[[1]], ms[[3]])), "one refine")
  V <- fe_space(ms[[2]])
  A <- assemble_matrix(V, "stiffness") + assemble_matrix(V, "mass")
  B <- cbind(assemble_vector(V, 1), assemble_vector(V, function(x, y) x))
  d <- solve_linear_system(A, B)
  expect_equal(dim(d$x), dim(B))
  expect_equal(d$info$status, "solved")
  for (pc in list("ssor", "none", function(r) 2 * r)) {
    s <- solve_linear_system(A, B[, 1], method = "cg", preconditioner = pc, tol = 1e-12)
    expect_lt(max(abs(s$x - d$x[, 1])), 1e-9)
  }
  s <- solve_linear_system(A, B, method = "gmres", preconditioner = "ssor", tol = 1e-12)
  expect_lt(max(abs(s$x - d$x)), 1e-9)
  expect_error(solve_linear_system(A, B[, 1], method = "cg", preconditioner = function(r) NA), "invalid")
})
