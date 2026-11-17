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
  m <- pde_mesh_rectangle(nx = 4)
  pf <- function(th) poisson_problem(m, diffusion = th[1], source = th[2], boundary = dirichlet_bc(0))
  obs <- pde_observations(c(.3, .5, .7), c(.5, .5, .3), rep(.1, 3))
  e <- estimate_parameters(pf, c(1, 1), obs, lower = .1, upper = 10, regularization = .01)
  id <- e$identifiability
  expect_identical(id$free_parameters, 2L)
  expect_lte(id$rank, 2L)
  expect_identical(ncol(id$null_space), 2L - id$rank)
  expect_length(id$singular_values, 2L)
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
