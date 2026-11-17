# pdecompR 0.1.0

First release. pdecompR is a finite-element framework for partial
differential equations on two-dimensional domains, written in R on top of
the sparse matrix classes of 'Matrix'.

## Geometry and meshes

* Rectangular and simple polygonal domains with validation.
* Structured triangular meshes and polygon meshes, with checks of
  connectivity, orientation and degeneracy, quality measures, named
  boundary markers and named regions (subdomains).
* Uniform refinement and newest vertex bisection; boundary markers and
  regions are inherited by refined meshes.
* Import of Gmsh 2.2 and 4.1 meshes (physical lines become boundary
  markers, physical surfaces become regions) and export to VTK, including
  quadratic P2 cells, displacement vectors and stress fields.

## Discretization

* Continuous Lagrange P1 and P2 elements on triangles, scalar and vector
  valued, with Gaussian quadrature on intervals and triangles.
* Sparse assembly of stiffness, mass, advection and boundary operators and
  load vectors. Terms with coefficients that are constant on each element
  are assembled from reference-element matrices, exactly and without a
  quadrature loop.
* Coefficients given as numbers, functions of `(x, y)` or `(x, y, t)`,
  finite-element functions, constant symmetric tensors or piecewise values
  by region (`region_coefficient()`).
* Dirichlet, Neumann and Robin conditions; pure Neumann problems with a
  compatibility check and a zero-mean gauge on each connected component.
* A weak-form language (`trial_function()`, `test_function()`, `grad()`,
  `dot()`, `integrate_domain()`, `integrate_boundary()`) for linear scalar
  problems.

## Problem classes

* Elliptic problems with diffusion, advection and reaction, with optional
  SUPG stabilization.
* Heat and general linear parabolic problems with backward Euler and
  Crank-Nicolson (Rannacher startup), time-dependent data and coefficients,
  and mass diagnostics.
* Reaction-diffusion systems with IMEX and fully implicit Newton time
  stepping.
* Quasilinear diffusion, semilinear reaction and p-Laplace problems solved
  by damped Newton iteration.
* Linear elasticity in plane strain and plane stress with displacement,
  roller and traction conditions, and strain, stress and von Mises output.

## Solvers

* Sparse Cholesky and LU with symmetric equilibration, and preconditioned
  conjugate gradients, BiCGSTAB and GMRES, with residual and conditioning
  diagnostics and an explicit status for every solve. Several right-hand
  sides are solved with one factorization, and time stepping warm starts
  iterative solvers.
* Jacobi, symmetric Gauss-Seidel (SSOR) and user-supplied preconditioners,
  and a geometric multigrid V-cycle on nested meshes
  (`multigrid_preconditioner()`) with mesh-independent iteration counts.

## Verification and adaptivity

* L2, H1 and energy error norms against exact or finite-element reference
  solutions on nested or unrelated meshes (integrated on the mesh
  overlay), and convergence studies with observed orders.
* Symbolic manufactured solutions through the optional 'symcompR' package.
* A residual a posteriori error estimator, Doerfler marking and an adaptive
  solve loop.
* A benchmark suite of ten reference problems in `benchmarks/`.

## Inference and control

* Discrete forward and adjoint sensitivities of linear elliptic problems,
  with finite-difference or exact operator derivatives
  (`affine_derivative()`) and Taylor remainder checks (`taylor_test()`).
* Parameter estimation with adjoint gradients, bounds (including fixed
  parameters), Gauss-Newton standard errors and identifiability
  diagnostics.
* Distributed optimal control with optional control bounds.
