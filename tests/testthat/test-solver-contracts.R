# Solver, constraint and periodic contracts.

test_that('incompatible coarse advection Neumann data is rejected', {
  p <- elliptic_problem(pde_mesh_rectangle(nx=1), advection=c(1,0), source=1)
  expect_error(solve_pde(p), 'incompatible', ignore.case=TRUE)
})

test_that('periodic eigenmodes enforce their constraints or explicitly reject support', {
  m <- pde_mesh_rectangle(nx=8)
  p <- elliptic_problem(m, degree=2, bcs=list(periodic_bc('left','right',c(1,0)),
                     dirichlet_bc(0,c('top','bottom'))))
  ev <- tryCatch(solve_eigen(p,n=4),error=identity)
  if(inherits(ev,'error')) {
    expect_match(conditionMessage(ev),'periodic',ignore.case=TRUE)
  } else {
    xy <- ev$vectors[[1]]$space$dof_coords
    l <- which(xy[,1]==0); r <- which(xy[,1]==1)
    l<-l[order(xy[l,2])];r<-r[order(xy[r,2])]
    expect_lt(max(vapply(ev$vectors,function(v)max(abs(v$values[l]-v$values[r])),numeric(1))),1e-8)
    expect_equal(ev$values/pi^2,c(1,4,5,5),tolerance=.01)
  }
})

test_that('periodic pairing respects scaling and translation', {
  for(cfg in list(c(1,0),c(1e-10,0),c(1,1e10))) {
    L<-cfg[1];origin<-cfg[2]
    m<-pde_mesh_rectangle(c(origin,origin+L),c(origin,origin+L),nx=4)
    p<-elliptic_problem(m,diffusion=L^2,reaction=1,degree=2,
       source=function(x,y)(pi^2+1)*sin(pi*(y-origin)/L),
       bcs=list(periodic_bc('left','right',c(L,0)),dirichlet_bc(0,c('top','bottom'))))
    s<-solve_pde(p);V<-s$space
    l<-boundary_dofs(V,'left');r<-boundary_dofs(V,'right')
    l<-l[order(V$dof_coords[l,2])];r<-r[order(V$dof_coords[r,2])]
    expect_lt(max(abs(s$u$values[l]-s$u$values[r])),1e-8)
    expect_lt(as.numeric(errornorm(s,function(x,y)sin(pi*(y-origin)/L)))/L,.005)
  }
})

test_that('periodic classes reject inconsistent prescribed values', {
  m<-pde_mesh_rectangle(nx=4)
  p<-elliptic_problem(m,bcs=list(periodic_bc('left','right',c(1,0)),
       dirichlet_bc(0,'left'),dirichlet_bc(1,'right')))
  expect_error(solve_pde(p),'inconsisten|conflict',ignore.case=TRUE)
})

test_that('reversed periodic duplicates preserve constraints or reject explicitly', {
  m<-pde_mesh_rectangle(nx=4)
  p<-elliptic_problem(m,reaction=1,source=function(x,y)x,bcs=list(
       periodic_bc('left','right',c(1,0)),periodic_bc('right','left',c(-1,0))))
  d<-tryCatch(discretize(p),error=identity)
  if(inherits(d,'error')) expect_match(conditionMessage(d),'periodic|cycle',ignore.case=TRUE) else {
    expect_equal(ncol(d$periodic),20L)
    s<-solve_pde(p);V<-s$space
    l<-boundary_dofs(V,'left');r<-boundary_dofs(V,'right')
    l<-l[order(V$dof_coords[l,2])];r<-r[order(V$dof_coords[r,2])]
    expect_equal(s$u$values[l],s$u$values[r],tolerance=1e-10)
  }
})

test_that('AMG survives uniform scaling of a positive definite system', {
  V<-fe_space(pde_mesh_rectangle(nx=10),2)
  A<-assemble_matrix(V,'stiffness')+assemble_matrix(V,'mass')
  b<-as.numeric(A%*%rep(1,nrow(A)))
  s<-solve_linear_system(1e-16*A,1e-16*b,method='cg',
                         preconditioner=amg_preconditioner(coarse_size=30))
  expect_identical(s$info$status,'converged')
  expect_lt(max(abs(s$x-1)),1e-8)
})

test_that('geometric multigrid permits an all-prescribed coarse mesh', {
  m<-pde_mesh_rectangle(nx=1);mf<-refine_uniform(m)
  p<-poisson_problem(mf,source=1,boundary=dirichlet_bc(0))
  s<-solve_pde(p,method='cg',preconditioner=multigrid_preconditioner(list(m,mf)))
  ref<-solve_pde(p)
  expect_identical(s$diagnostics$solver$status,'converged')
  expect_equal(s$u$values,ref$u$values,tolerance=1e-10)
})

test_that('normal-only closed Stokes constraints receive a pressure gauge', {
  m<-pde_mesh_rectangle(nx=4)
  bs<-list(velocity_bc(0,c('left','right'),components='x'),
           velocity_bc(0,c('bottom','top'),components='y'))
  s<-solve_pde(stokes_problem(m,body_force=c(1,0),bcs=bs))
  expect_lt(max(abs(c(s$velocity$ux$values,s$velocity$uy$values))),1e-9)
  expect_lt(as.numeric(errornorm(s$pressure,function(x,y)x-.5)),1e-9)
  expect_identical(s$diagnostics$pressure_gauge,'zero mean')
})

test_that("the automatic solver choice does not depend on the scale of A", {
  A <- matrix(c(2, 1, 0, 2), 2); b <- c(2, 3)
  for (k in c(1, 1e-15, 1e15)) {
    s <- solve_linear_system(k * A, k * b)
    expect_equal(s$x, c(1, 1), tolerance = 1e-12)
    expect_equal(s$info$method, "lu")
  }
  # A small asymmetric block next to a large symmetric one is still asymmetric.
  A <- diag(c(1, 1e-13, 1e-13)); A[3, 2] <- 0.9e-13
  s <- solve_linear_system(A, c(0, 1e-13, 1e-13))
  expect_equal(s$x, c(0, 1, 0.1), tolerance = 1e-12)
})

test_that("Matrix right-hand sides, flexible GMRES and AMG scaling", {
  s <- solve_linear_system(diag(c(2, 3)), Matrix::Matrix(matrix(1:4, 2)))
  expect_equal(s$x, cbind(c(.5, 1 / 3) * 1:2, c(1.5, 4 / 3)), tolerance = 1e-12)
  # A nonlinear preconditioner: flexible GMRES applies it consistently.
  s <- solve_linear_system(matrix(1), 1, method = "gmres", maxit = 1, preconditioner = function(r) r * (1 + r^2))
  expect_equal(s$x, 1, tolerance = 1e-12)
  m <- pde_mesh_rectangle(nx = 24)
  it <- function(k) solve_pde(poisson_problem(m, diffusion = k, source = k, bcs = dirichlet_bc(0)), method = "cg",
                              preconditioner = amg_preconditioner(coarse_size = 50))$diagnostics$solver$iterations
  expect_equal(it(1e160), it(1))
})

test_that("input contracts: tensors, Robin, eigenpairs, polygons, Gmsh, bounds", {
  m <- pde_mesh_rectangle(nx = 2)
  s <- solve_pde(poisson_problem(m, diffusion = tensor_coefficient(1e-180), source = 1e-180, boundary = dirichlet_bc(0)))
  expect_equal(s$u$values, solve_pde(poisson_problem(m, source = 1, boundary = dirichlet_bc(0)))$u$values, tolerance = 1e-12)
  expect_error(solve_pde(poisson_problem(m, boundary = robin_bc(function(x, y) -1, 1))), "non-negative")
  e <- solve_eigen(poisson_problem(m, boundary = dirichlet_bc(0)), n = 1)
  expect_length(e$values, 1)
  V <- fe_space(pde_mesh_rectangle(nx = 3)); N <- V$ndof - length(boundary_dofs(V))
  expect_length(solve_eigen(poisson_problem(pde_mesh_rectangle(nx = 3), boundary = dirichlet_bc(0)), n = N)$values, N)
  expect_equal(integrate_mesh(pde_mesh_polygon(1e15 + rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))), 1), 1)
  f <- tempfile(fileext = ".msh")
  writeLines(c("$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$Nodes", "3", "1 0 0 0", "2 1 0 0", "3 0 1 0", "$EndNodes",
               "$Elements", "2", "1 1 2 7 1 1 99", "2 2 0 1 2 3", "$EndElements"), f)
  expect_error(read_gmsh(f), "line element references a node")
  expect_error(optimal_control(m, target = 1, alpha = .01, lower = 1, upper = -1), "must not exceed")
})

test_that("fixed parameters keep their regularization in the objective", {
  m <- pde_mesh_rectangle(nx = 2)
  pf <- function(z) poisson_problem(m, source = z[2], boundary = dirichlet_bc(0))
  e <- estimate_parameters(pf, c(2, 0), pde_observations(.5, .5, 0), lower = c(2, -1), upper = c(2, 1),
                           regularization = 1, theta_ref = c(0, 0))
  expect_equal(e$objective, 2, tolerance = 1e-8)
})

test_that("gauges are decided per component after the constraints", {
  P <- rbind(c(0, 0), c(1, 0), c(0, 1), c(2, 0), c(3, 0), c(2, 1))
  b <- set_boundary_marker(pde_mesh(P, rbind(1:3, 4:6)), "fixed", function(x, y) x < 1.5)
  s <- solve_pde(poisson_problem(b, boundary = dirichlet_bc(1, "fixed")))
  expect_equal(s$u$values, c(1, 1, 1, 0, 0, 0), tolerance = 1e-12)
  m2 <- pde_mesh_rectangle(nx = 2)
  expect_error(solve_pde(poisson_problem(m2, source = 1e-310)), "incompatible")
  expect_error(solve_pde(poisson_problem(m2, source = 1), compatibility_tol = .02), "at most 0.01")
  s <- solve_pde(poisson_problem(m2, source = function(x, y) x - .5 + 1e-9))
  expect_false(s$diagnostics$gauge$original_problem_solved)
  bc <- list(periodic_bc("left", "right", c(1, 0)), dirichlet_bc(function(x, y) x, "bottom"), dirichlet_bc(1e12, "top"))
  expect_error(solve_pde(poisson_problem(m2, boundary = bc)), "Inconsistent prescribed values")
  pp <- function(a) poisson_problem(m2, bcs = list(periodic_bc("left", "right", c(1, 0)), dirichlet_bc(a, c("bottom", "top"))))
  expect_equal(range(pde_sensitivity(pp, 1, derivative = affine_derivative(pp, 1))$sensitivities), c(1, 1), tolerance = 1e-10)
})

test_that("Stokes reports free translations and too few velocity unknowns", {
  m <- pde_mesh_rectangle(nx = 2)
  expect_error(solve_pde(stokes_problem(m, body_force = c(1, 0))), "not constrained")
  expect_error(solve_pde(stokes_problem(pde_mesh_rectangle(nx = 1), bcs = velocity_bc(0))), "pressure is not determined")
  n <- nrow(m$nodes); e <- cbind(m$boundary$n1, m$boundary$n2)
  Q <- m$nodes * 2^-30; Q[, 1] <- 2 + Q[, 1]
  b <- pde_mesh(rbind(m$nodes, Q), rbind(m$triangles, m$triangles + n), rbind(e, e + n),
                c(rep("wall", nrow(e)), ifelse(m$boundary$marker == "right", "out", "wall")))
  s <- suppressWarnings(solve_pde(stokes_problem(b, bcs = list(velocity_bc(0, "wall"), traction_bc(c(-1, 0), "out")))))
  expect_equal(range(tail(s$pressure$values, n)), c(1, 1), tolerance = 1e-8)
})

test_that("repeated assembly switches to plans without changing results", {
  V <- fe_space(pde_mesh_rectangle(nx = 4), 2)
  geo <- element_geometry(V, NULL)
  loc <- local_bilinear(geo, diff = coef_qp(function(x, y) 1 + x, geo))
  lv <- local_linear(geo, f = coef_qp(function(x, y) y, geo))
  asm <- repeated_assembler(V)
  for (i in 1:25) { A <- asm$matrix(loc); b <- asm$vector(lv) }
  expect_lt(max(abs(A - global_matrix(V, loc))), 1e-14)
  expect_lt(max(abs(b - global_vector(V, lv))), 1e-15)
  # The plan closure keeps only its summation map and pattern.
  expect_setequal(ls(environment(assembly_plan(V))), c("S", "tmpl"))
})

test_that("errors against a reference on the same space use the difference directly", {
  s <- solve_pde(poisson_problem(pde_mesh_rectangle(nx = 6), degree = 2, source = 1, boundary = dirichlet_bc(0)))
  v <- interpolate(function(x, y) x * y, s$space); z <- s$u$values - v$values
  expect_equal(as.numeric(errornorm(s, v)), sqrt(sum(z * as.numeric(assemble_matrix(s$space, "mass") %*% z))),
               tolerance = 1e-12)
  expect_equal(as.numeric(errornorm(s, v, "H1_semi")),
               sqrt(sum(z * as.numeric(assemble_matrix(s$space, "stiffness") %*% z))), tolerance = 1e-12)
})

test_that("moving Dirichlet data do not reassemble a static load", {
  calls <- 0L
  f <- function(x, y) { calls <<- calls + 1L; x }
  solve_pde(heat_problem(pde_mesh_rectangle(nx = 3), source = f, bcs = dirichlet_bc(function(x, y, t) t * x),
                         time = c(0, .1)), dt = .01)
  expect_equal(calls, 1L)
})
