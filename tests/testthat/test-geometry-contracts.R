# Geometry, integration and verification contracts.

test_that("overlay error retains real fine elements and is symmetric", {
  coarse <- pde_mesh(rbind(c(0,0),c(1,0),c(1,1),c(0,1),c(.37,.57)),
                     rbind(c(1,2,5),c(2,3,5),c(3,4,5),c(4,1,5)))
  fine <- pde_mesh_rectangle(nx=2,pattern="left")
  center <- which(fine$nodes[,1]==.5 & fine$nodes[,2]==.5)
  for(i in 1:20) fine <- refine_mesh(fine,which(rowSums(fine$triangles==center)>0))
  V <- fe_space(fine); q <- numeric(V$ndof);q[center]<-1
  M <- assemble_matrix(V,"mass")
  q <- q/sqrt(as.numeric(crossprod(q,M%*%q)))
  u <- interpolate(0,fe_space(coarse));v <- fe_function(V,q)
  expect_gte(min(mesh_quality(fine)$min_angle),44.99)
  expect_equal(as.numeric(errornorm(u,v)),1,tolerance=1e-8)
  expect_equal(as.numeric(errornorm(u,v)),as.numeric(errornorm(v,u)),tolerance=1e-8)
  bounded <- numeric(V$ndof);bounded[center] <- 1;hat <- fe_function(V,bounded)
  K <- assemble_matrix(V,"stiffness")
  expect_equal(as.numeric(errornorm(u,hat,"H1_semi")),
               sqrt(as.numeric(crossprod(bounded,K%*%bounded))),tolerance=1e-8)
})

test_that("default verification gradient follows domain scale after translation", {
  z <- 1e6
  V <- fe_space(pde_mesh_rectangle(c(z,z+1),c(0,1),nx=16),2)
  f <- function(x,y) sin(x-z);g <- function(x,y)cbind(cos(x-z),0*y)
  u <- interpolate(f,V)
  oracle <- as.numeric(errornorm(u,f,"H1_semi",gradient=g))
  default <- as.numeric(errornorm(u,f,"H1_semi"))
  expect_lt(oracle,2e-4)
  expect_lt(abs(default/oracle-1),.02)
})

test_that("contained but different domains do not bypass integration contract", {
  big <- fe_space(pde_mesh_rectangle(nx=1))
  small <- fe_space(pde_mesh_rectangle(c(0,.5),c(0,.5),nx=1))
  a <- interpolate(0,small);b <- interpolate(1,big)
  expect_error(errornorm(a,b),"domain")
  expect_error(errornorm(b,a),"domain")
  expect_error(l2_project(b,small),"domain")
})

test_that("curved hole refinement rejects common-domain comparisons", {
  th <- 2*pi*(1:8)/8
  d <- pde_domain("polygon",vertices=2*cbind(cos(th),sin(th)),markers="outer",
    holes=list(cbind(cos(th),sin(th))),hole_markers="inner")
  m <- curved_boundary(pde_mesh_polygon(d),"inner",function(x,y)cbind(x,y)/sqrt(x*x+y*y))
  fine <- refine_uniform(m)
  a <- interpolate(1,fe_space(m));b <- interpolate(0,fe_space(fine))
  expect_gt(abs(sum(mesh_quality(m)$area)-sum(mesh_quality(fine)$area)),.2)
  expect_error(errornorm(a,b),"domain")
  expect_error(l2_project(a,b$space),"domain")
})

test_that("strict overlap validation detects a contained disconnected component", {
  P <- rbind(c(0,0),c(4,0),c(0,4),c(.5,.5),c(1,.5),c(.5,1))
  expect_error(pde_mesh(P,rbind(c(1,2,3),c(4,5,6))),"overlap")
})

test_that("tensor energy integration does not evaluate outside physical domain", {
  V <- fe_space(pde_mesh_rectangle(c(2,3),c(2,3),nx=2))
  u <- interpolate(function(x,y)x,V)
  K <- tensor_coefficient(function(x,y){
    if(any(x<2|x>3|y<2|y>3))stop("outside domain")
    rep(1,length(x))
  })
  expect_true(all(is.finite(assemble_matrix(V,"stiffness",coef=K)@x)))
  expect_equal(as.numeric(errornorm(u,function(x,y)0,"energy",
    gradient=function(x,y)cbind(0*x,0*y),diffusion=K)),1,tolerance=1e-12)
})

test_that("hole containment uses local geometry rather than absolute origin", {
  o <- rbind(c(0,-2),c(2,0),c(0,2),c(-2,0));h <- .25*o
  expect_no_error(pde_domain("polygon",vertices=1e-7*o,holes=list(1e-7*h)))
  expect_no_error(pde_domain("polygon",vertices=5e5+1e-4*o,holes=list(5e5+1e-4*h)))
})

test_that("relabeling generic boundary drops its curve as documented", {
  th <- 2*pi*(1:8)/8
  m <- curved_boundary(pde_mesh_polygon(cbind(cos(th),sin(th))),"boundary",
                       function(x,y)cbind(x,y)/sqrt(x*x+y*y))
  m <- set_boundary_marker(m,"straight",function(x,y)rep(TRUE,length(x)))
  fine <- refine_uniform(m);n <- nrow(m$nodes);E <- m$edges
  midpoint <- (m$nodes[E[,1],]+m$nodes[E[,2],])/2
  expect_equal(fine$nodes[n+seq_len(nrow(E)),],midpoint,tolerance=1e-12)
})

test_that("curves follow their edges through local refinement and partial relabeling", {
  th <- 2 * pi * (1:8) / 8
  m <- curved_boundary(pde_mesh_polygon(cbind(cos(th), sin(th))), "boundary",
                       function(x, y) cbind(x, y) / sqrt(x * x + y * y))
  near <- function(f) which(rowSums(matrix(sqrt(rowSums(f$nodes[f$triangles, ]^2)) > 0.99, ncol = 3)) >= 2)
  f <- refine_mesh(refine_mesh(m, near(m)[1:2]), near(m)[1:3])
  bd <- f$boundary
  expect_false(anyNA(bd$curve))
  expect_equal(sqrt(rowSums(f$nodes[unique(c(bd$n1, bd$n2)), ]^2)), rep(1, length(unique(c(bd$n1, bd$n2)))),
               tolerance = 1e-12)
  # Relabel the edges in the upper half straight; only the lower half stays curved.
  s <- set_boundary_marker(m, "flat", function(x, y) y > 1e-9)
  r <- refine_uniform(s); n <- nrow(s$nodes); E <- s$edges
  rad <- sqrt(rowSums(r$nodes[n + seq_len(nrow(E)), ]^2))
  bmid <- match(paste(pmin(s$boundary$n1, s$boundary$n2), pmax(s$boundary$n1, s$boundary$n2)),
                paste(pmin(E[, 1], E[, 2]), pmax(E[, 1], E[, 2])))
  upper <- s$boundary$marker == "flat"
  expect_equal(rad[bmid[!upper]], rep(1, sum(!upper)), tolerance = 1e-12)
  expect_true(all(rad[bmid[upper]] < 1 - 1e-3))
})

test_that("nesting is certified per parent and thin crossing strips are integrated", {
  h <- 2^-40
  kink <- function(a) {
    b <- pde_mesh_rectangle(nx = 2, ny = 1); P <- b$nodes; P[P[, 1] == .5, 1] <- a
    interpolate(function(x, y) pmax(0, x - a) / sqrt(h), fe_space(pde_mesh(P, b$triangles)))
  }
  u <- kink(.5); v <- kink(.5 + h)
  expect_equal(as.numeric(errornorm(u, v, "H1_semi")), 1, tolerance = 1e-6)
  expect_equal(as.numeric(errornorm(v, u, "H1_semi")), 1, tolerance = 1e-6)
})

test_that("a missing small component is a different domain", {
  m <- pde_mesh_rectangle(nx = 2); h <- 2^-16; n <- nrow(m$nodes)
  b <- pde_mesh(rbind(m$nodes, c(2, 0), c(2 + h, 0), c(2, h)), rbind(m$triangles, n + 1:3))
  u <- fe_function(fe_space(b), c(rep(0, n), rep(2 / h^2, 3)))
  z <- interpolate(0, fe_space(m))
  expect_error(l2_project(u, z$space), "do not cover the same domain")
  expect_error(errornorm(u, z), "do not cover the same domain")
  expect_error(errornorm(z, u), "do not cover the same domain")
})

test_that("components touching another element's interior edge overlap", {
  sq <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  P <- rbind(sq, c(.375, .375), c(.625, .375), c(.5, .75))
  expect_error(pde_mesh(P, rbind(c(1, 2, 3), c(1, 3, 4), 5:7)), "overlap")
})

test_that("curved refinement validates the projected mesh", {
  D <- rbind(c(1, 0), c(0, 1), c(-1, 0), c(0, -1)); T <- rbind(c(1, 2, 3), c(1, 3, 4))
  b <- pde_mesh(rbind(D, D + 1.25), rbind(T, T + 4))
  b <- set_boundary_marker(b, "one", function(x, y) x + y < 1.2)
  b <- set_boundary_marker(b, "two", function(x, y) x + y > 1.2)
  circle <- function(z) function(x, y) z + cbind(x - z, y - z) / sqrt((x - z)^2 + (y - z)^2)
  b <- curved_boundary(curved_boundary(b, "one", circle(0)), "two", circle(1.25))
  expect_error(refine_uniform(b), "invalid mesh")
  expect_error(refine_mesh(b, 1:4), "invalid mesh")
})

test_that("overlay projection and estimators are translation invariant", {
  z <- function(o) {
    V <- function(n, p) fe_space(pde_mesh_rectangle(x = o + c(0, 1), y = o + c(0, 1), nx = n, pattern = p))
    a <- V(3, "right"); b <- V(4, "left")
    f <- interpolate(function(x, y) x - o, a)
    max(abs(l2_project(f, b)$values - (b$dof_coords[, 1] - o)))
  }
  expect_lt(max(sapply(c(0, 1e9, 1e13), z)), 1e-12)
  for (o in c(0, 2^36)) {
    b <- pde_mesh_rectangle(x = c(o, o + 1), nx = 4)
    s <- solve_pde(poisson_problem(b, boundary = dirichlet_bc(function(x, y) x - o)))
    expect_lt(estimate_error(s)$total, 1e-10)
  }
  e <- function(o) {
    m <- set_region(pde_mesh_rectangle(x = o + c(0, 1), nx = 2, ny = 1), "high", function(x, y) x - o > .5)
    k <- region_coefficient(m, c(domain = 1, high = 2))
    estimate_error(solve_pde(poisson_problem(m, diffusion = k, bcs = list(dirichlet_bc(0, "left"),
                                                                         neumann_bc(1, "right")))))$total
  }
  expect_lt(e(1e10), 1e-10)
})

test_that("polygons with several holes are meshed", {
  sq <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  holes <- lapply(c(1, 3), function(y) sweep(sq, 2, c(3, y), "+"))
  m <- pde_mesh_polygon(pde_domain("polygon", vertices = 10 * sq, holes = holes))
  expect_equal(integrate_mesh(m, 1), 98, tolerance = 1e-12)
  holes3 <- lapply(c(1, 4, 7), function(y) sweep(sq, 2, c(y, y), "+"))
  m3 <- pde_mesh_polygon(pde_domain("polygon", vertices = 10 * sq, holes = holes3))
  expect_equal(integrate_mesh(m3, 1), 97, tolerance = 1e-12)
})
