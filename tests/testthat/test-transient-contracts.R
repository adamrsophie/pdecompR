# Derivative, restart and time contracts.

test_that("transient SUPG adjoints include parameter-dependent mass derivatives", {
  m <- pde_mesh_rectangle(nx=8)
  ob <- function(u) list(value=sum(u^2)/2,gradient=u)
  pf <- function(th) parabolic_problem(m,diffusion=.01,advection=c(th,.2),source=0,
    initial=function(x,y)sin(pi*x)*sin(pi*y),bcs=dirichlet_bc(0),time=c(0,.2),stabilization="supg")
  for(me in c("backward-euler","crank-nicolson")) {
    J <- function(th) ob(solve_pde(pf(th),dt=.02,method=me,startup=0)$state$u)$value
    fd <- (J(1+1e-5)-J(1-1e-5))/2e-5
    g <- transient_sensitivity(pf,1,ob,dt=.02,method=me)$gradient
    expect_equal(g,fd,tolerance=1e-6,info=me)
  }
})

test_that("periodic sensitivities use the same constrained state as solve_pde", {
  m <- pde_mesh_rectangle(nx=8)
  bc <- list(periodic_bc("left","right",c(1,0)),dirichlet_bc(0,c("bottom","top")))
  pf <- function(th) elliptic_problem(m,source=function(x,y)th*sin(2*pi*x),bcs=bc)
  ob <- function(u)list(value=sum(u^2)/2,gradient=u)
  ss <- pde_sensitivity(pf,1,objective=ob,method="adjoint")
  exact <- solve_pde(pf(1))
  expect_equal(ss$solution$u$values,exact$u$values,tolerance=1e-10)
  fd <- (ob(solve_pde(pf(1+1e-5))$u$values)$value-ob(solve_pde(pf(1-1e-5))$u$values)$value)/2e-5
  expect_equal(ss$gradient,fd,tolerance=1e-6)
})

test_that("periodic optimal control satisfies its own state equation", {
  m <- pde_mesh_rectangle(nx=8)
  bc <- list(periodic_bc("left","right",c(1,0)),dirichlet_bc(0,c("bottom","top")))
  r <- optimal_control(m,target=function(x,y)sin(2*pi*x)*sin(pi*y),alpha=.001,bcs=bc)
  s <- solve_pde(elliptic_problem(m,source=r$control,bcs=bc))
  expect_equal(r$state$values,s$u$values,tolerance=1e-8)
})

test_that("restart checks or maps the finite-element numbering", {
  m <- pde_mesh_rectangle(nx=4)
  perm <- rev(seq_len(nrow(m$nodes))); inv <- match(seq_len(nrow(m$nodes)),perm)
  m2 <- pde_mesh(m$nodes[perm,,drop=FALSE],matrix(inv[m$triangles],ncol=3))
  pf <- function(mm,T)heat_problem(mm,initial=function(x,y)x,bcs=dirichlet_bc(function(x,y)x),time=c(0,T))
  a <- solve_pde(pf(m,.01),dt=.001,method="backward-euler")
  r <- tryCatch(solve_pde(pf(m2,.0101),dt=.0001,method="backward-euler",restart=a),error=identity)
  if(inherits(r,"error")) expect_match(conditionMessage(r),"restart|mesh|space|numbering")
  else expect_equal(r$state$u,r$space$dof_coords[,1],tolerance=1e-10)
})

test_that("restart checks or maps reaction-diffusion species identities", {
  m <- pde_mesh_rectangle(nx=2)
  pf <- function(sp,T)reaction_diffusion_problem(m,species=sp,diffusion=0,
    reaction=function(u,x,y,t)0*u,initial=setNames(as.list(seq_along(sp)),sp),time=c(0,T))
  a <- solve_pde(pf(c("A","B"),.1),dt=.1)
  r <- tryCatch(solve_pde(pf(c("B","A"),.2),dt=.1,restart=a),error=identity)
  if(inherits(r,"error")) expect_match(conditionMessage(r),"restart|species")
  else {
    expect_equal(unname(r$state$u[,"A"]),rep(1,nrow(r$state$u)),tolerance=1e-10)
    expect_equal(unname(r$state$u[,"B"]),rep(2,nrow(r$state$u)),tolerance=1e-10)
  }
})

test_that("identifiability detects fewer observations than free parameters", {
  m <- pde_mesh_rectangle(nx=4)
  pf <- function(th)poisson_problem(m,diffusion=th[1],source=th[2],boundary=dirichlet_bc(0))
  e <- estimate_parameters(pf,c(1,1),pde_observations(.5,.5,.1),lower=.1,upper=10,regularization=.01)
  expect_true(is.infinite(e$identifiability$condition))
})

test_that("a constant reaction does not falsely certify Neumann uniqueness", {
  p <- nonlinear_problem(pde_mesh_rectangle(nx=3),reaction=function(u)rep(1,length(u)),
    reaction_derivative=function(u)0*u,source=1,initial=7)
  r <- tryCatch(solve_pde(p),error=identity)
  if(inherits(r,"error")) expect_match(conditionMessage(r),"unique|Neumann|nullspace|gauge")
  else expect_false(is.null(r$diagnostics$gauge))
})

test_that("the zero-Jacobian cubic Neumann solution is not falsely rejected", {
  p <- nonlinear_problem(pde_mesh_rectangle(nx=3),reaction=function(u)u^3,
    reaction_derivative=function(u)3*u^2,source=0,initial=0)
  r <- solve_pde(p)
  expect_true(r$converged)
  expect_equal(r$u$values,rep(0,r$space$ndof),tolerance=1e-12)
})

test_that("adaptive time integration reaches its endpoint at small time units", {
  T <- 1e-13
  p <- heat_problem(pde_mesh_rectangle(nx=2),source=1/T,initial=0,time=c(0,T))
  r <- solve_pde(p,dt=T/10,method="sdirk2",adaptive=TRUE,tol=1e-6)
  expect_gt(r$diagnostics$steps,0)
  expect_equal(r$state$u,rep(1,r$space$ndof),tolerance=1e-8)
  expect_equal(tail(r$times,1),T,tolerance=1e-25)
})

test_that("one-step waves do not query boundary data after the end time", {
  p <- wave_problem(pde_mesh_rectangle(nx=2),source=2,initial=0,velocity=0,
    bcs=dirichlet_bc(function(x,y,t){if(t>1)stop("outside requested time interval");t^2}),time=c(0,1))
  expect_no_error(solve_pde(p,dt=1))
})

test_that("identifiability is full rank with enough informative observations", {
  m <- pde_mesh_rectangle(nx = 6)
  pf <- function(th) elliptic_problem(m, diffusion = th[1], reaction = th[2], source = 1, bcs = dirichlet_bc(0))
  xo <- c(.2, .5, .7, .4); yo <- c(.5, .5, .3, .8)
  truth <- solve_pde(pf(c(1, 5)))
  obs <- pde_observations(xo, yo, evaluate(truth$u, xo, yo), sd = 0.01)
  e <- estimate_parameters(pf, c(1.5, 3), obs, lower = c(.1, .1), upper = c(10, 50))
  id <- e$identifiability
  expect_identical(id$rank, 2L)
  expect_true(is.finite(id$condition))
  expect_identical(ncol(id$null_space), 0L)
  expect_false(is.null(id$correlation))
})

test_that("transient SUPG adjoint gradients hold for diffusion and with checkpoints", {
  m <- pde_mesh_rectangle(nx = 8)
  ob <- function(u) list(value = sum(u^2) / 2, gradient = u)
  pf <- function(th) parabolic_problem(m, diffusion = .01 * th[2], advection = c(th[1], .2), source = 0,
    initial = function(x, y) sin(pi * x) * sin(pi * y), bcs = dirichlet_bc(0), time = c(0, .2),
    stabilization = "supg")
  for (me in c("backward-euler", "crank-nicolson")) {
    J <- function(th) ob(solve_pde(pf(th), dt = .02, method = me, startup = 0)$state$u)$value
    fd <- sapply(1:2, function(j) { e <- c(0, 0); e[j] <- 1e-5; (J(c(1, 1) + e) - J(c(1, 1) - e)) / 2e-5 })
    g <- transient_sensitivity(pf, c(1, 1), ob, dt = .02, method = me, checkpoint = 3L)$gradient
    expect_equal(g, fd, tolerance = 1e-6, info = me)
  }
})

test_that("old restart states without a layout version are rejected", {
  p <- heat_problem(pde_mesh_rectangle(nx = 2), initial = 0, time = c(0, .1))
  a <- solve_pde(p, dt = .05)
  a$state$version <- NULL
  expect_error(solve_pde(heat_problem(pde_mesh_rectangle(nx = 2), initial = 0, time = c(0, .2)), dt = .05,
                         restart = a), "state layout")
})

test_that("Newmark startup stays second order with moving boundary data", {
  m <- pde_mesh_rectangle(nx = 4); ex <- function(t) sin(2 * t) + t^3
  err <- sapply(c(.1, .05, .025), function(h) {
    p <- wave_problem(m, source = function(x, y, t) (-4 * sin(2 * t) + 6 * t) * x, initial = 0,
                      velocity = function(x, y) 2 * x, bcs = dirichlet_bc(function(x, y, t) ex(t) * x),
                      time = c(0, 1))
    s <- solve_pde(p, dt = h)
    max(abs(s$state$u - ex(1) * s$space$dof_coords[, 1]))
  })
  expect_gt(min(log2(err[-3] / err[-1])), 1.9)
})

test_that("time intervals too short to represent are rejected", {
  p <- heat_problem(pde_mesh_rectangle(nx = 2), initial = 0, time = c(1e6, 1e6 + 1e-9))
  expect_error(solve_pde(p, dt = 1e-10, method = "sdirk2", adaptive = TRUE), "too short")
})

test_that("time-dependent tensors are reassembled and adaptive SUPG uses the Galerkin norm", {
  m <- pde_mesh_rectangle(nx = 2)
  k <- function(x, y, t) 1 + 8 * t
  step <- function(d) solve_pde(heat_problem(m, diffusivity = d, time = c(0, .125), bcs = dirichlet_bc(0),
                                             initial = function(x, y) 16 * x * (1 - x) * y * (1 - y)),
                                dt = .125, method = "backward-euler")$state$u[5]
  expect_equal(step(tensor_coefficient(k)), step(k), tolerance = 1e-12)
  expect_equal(step(k), 1 / 9, tolerance = 1e-12)
  g <- 1 - 1 / sqrt(2); h <- .01
  p <- parabolic_problem(pde_mesh_rectangle(nx = 1), diffusion = .001, advection = c(1, 0), stabilization = "supg",
                         source = function(x, y, t) (t / h - g) * (1 - x), time = c(0, h))
  expect_gt(solve_pde(p, dt = h, method = "sdirk2", adaptive = TRUE, tol = 1e-8)$diagnostics$steps, 1)
})

test_that("restart checks connectivity and rebalances the acceleration", {
  pp <- function(pattern, T) heat_problem(pde_mesh_rectangle(nx = 1, pattern = pattern), initial = function(x, y) x * y,
                                          time = c(0, T))
  a <- solve_pde(pp("right", .1), dt = .1)
  expect_error(solve_pde(pp("left", .2), dt = .1, restart = a), "different finite-element space")
  m <- pde_mesh_rectangle(nx = 2)
  a <- solve_pde(wave_problem(m, source = 0, time = c(0, 1)), dt = 1)
  b <- solve_pde(wave_problem(m, source = 1, time = c(0, 2)), dt = 1, restart = a)
  expect_equal(c(b$state$u[5], b$state$v[5], b$state$a[5]), c(.5, 1, 1), tolerance = 1e-12)
  rd <- function(initial, T) reaction_diffusion_problem(m, initial = initial, reaction = function(u, x, y, t) 0 * u,
                                                        time = c(0, T))
  a <- solve_pde(rd(1, .5), dt = .5)
  r <- solve_pde(rd(function(x, y) stop("initial executed"), 1), dt = .5, restart = a)
  expect_equal(range(r$state$u), c(1, 1))
  q <- update_mesh(elastodynamics_problem(elasticity_problem(m, bcs = displacement_bc(0, "left"))), refine_uniform(m))
  expect_equal(solve_pde(q, dt = 1)$diagnostics$ndof, 2 * nrow(q$mesh$nodes))
})

test_that("time grids are built in elapsed time with exact endpoints", {
  m <- pde_mesh_rectangle(nx = 2); o <- 2^40
  p <- heat_problem(m, source = 1, time = c(o, o + 1))
  expect_equal(solve_pde(p, dt = .25)$times - o, c(0, .25, .5, .75, 1))
  expect_equal(solve_pde(p, dt = .25, method = "sdirk2", adaptive = TRUE)$state$t, o + 1)
  expect_error(solve_pde(heat_problem(m, source = 1, time = c(o, o + 2^-12)), dt = 2^-14), "too small")
  g <- function(x, y, t) if (t > .1) stop("past T") else 0
  w <- solve_pde(wave_problem(pde_mesh_rectangle(nx = 1), time = c(-1, .1), bcs = dirichlet_bc(g)), dt = 2)
  expect_identical(w$times[2], .1)
  expect_identical(w$state$t, .1)
})

test_that("the transient adjoint owns its trajectory and checks every time level", {
  m <- pde_mesh_rectangle(nx = 2)
  pf <- function(z) parabolic_problem(m, reaction = z, initial = 1)
  J <- function(u) list(value = mean(u), gradient = rep(1 / length(u), length(u)))
  g <- function(...) transient_sensitivity(pf, 1, J, dt = .25, ...)$gradient
  expect_equal(g(), -1.25^-5, tolerance = 1e-7)
  expect_equal(g(output_times = (1:4) / 4), -1.25^-5, tolerance = 1e-7)
  expect_equal(g(output_times = 1, checkpoint = 2), -1.25^-5, tolerance = 1e-7)
  pf2 <- function(z) heat_problem(m, bcs = dirichlet_bc(function(x, y, t) z * t * (t - .5) * (t - 1)))
  J5 <- function(u) list(value = u[5], gradient = as.numeric(seq_along(u) == 5))
  expect_error(transient_sensitivity(pf2, 1, J5, dt = .25), "must not depend on theta")
})

test_that("Newton convergence floors are dimensionless", {
  m <- pde_mesh_rectangle(nx = 2)
  for (k in c(1, 1e-18)) {
    s <- solve_pde(nonlinear_problem(m, diffusion = k, source = k, bcs = dirichlet_bc(0), initial = 0))
    expect_equal(s$u$values[5], .0625, tolerance = 1e-10)
  }
  for (L in c(1, 1e-6)) {
    b <- pde_mesh_rectangle(c(0, L), c(0, L), nx = 2)
    p <- reaction_diffusion_problem(b, diffusion = 0, reaction = function(u, x, y, t) -u, initial = 1)
    expect_equal(solve_pde(p, dt = 1, method = "implicit-euler")$state$u[5], .5, tolerance = 1e-10)
  }
})

test_that("nonlinear initialization and the flat constant mode", {
  m <- pde_mesh_rectangle(nx = 2)
  p <- nonlinear_problem(m, reaction = function(u) u^3, reaction_derivative = function(u) 3 * u^2, source = 1)
  expect_equal(range(solve_pde(p)$u$values), c(1, 1), tolerance = 1e-8)
  gf <- function(u) pmax(u - 2.25, 0)^3 - pmax(1.75 - u, 0)^3
  expect_error(solve_pde(nonlinear_problem(pde_mesh_rectangle(nx = 1), reaction = gf, initial = 2)), "not unique")
  p0 <- nonlinear_problem(m, reaction = function(u) u^3, reaction_derivative = function(u) 3 * u^2, source = 0)
  expect_equal(max(abs(solve_pde(p0)$u$values)), 0)
})

test_that("numerical gradients resolve translated coordinates and projections integrate the mass exactly", {
  z <- function(o) {
    V <- fe_space(pde_mesh_rectangle(x = o + c(0, 1), nx = 8))
    errornorm(interpolate(0, V), function(x, y) sin(2 * pi * (x - o)), "H1_semi")
  }
  expect_equal(z(1e12), z(0), tolerance = 1e-4)
  V <- fe_space(pde_mesh_rectangle(nx = 1), 2)
  expect_equal(l2_project(1, V, quad_degree = 2)$values, rep(1, V$ndof), tolerance = 1e-12)
})
