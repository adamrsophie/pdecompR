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
