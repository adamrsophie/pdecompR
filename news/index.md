# Changelog

## pdecompR 0.1.0

First release. pdecompR is a finite-element framework for partial
differential equations on two-dimensional domains, written in R on top
of the sparse matrix classes of ‘Matrix’.

### Geometry and meshes

- Rectangular and simple polygonal domains with validation.
- Structured triangular meshes and polygon meshes, with checks of
  connectivity, orientation and degeneracy, quality measures, named
  boundary markers and named regions (subdomains).
- Uniform refinement and newest vertex bisection; boundary markers and
  regions are inherited by refined meshes.
- Import of Gmsh 2.2 and 4.1 meshes (physical lines become boundary
  markers, physical surfaces become regions) and export to VTK,
  including quadratic P2 cells, displacement vectors and stress fields,
  in legacy and XML (`.vtu`) formats, and ParaView time series
  ([`write_pvd()`](https://adamrsophie.github.io/pdecompR/reference/write_pvd.md)).

### Discretization

- Continuous Lagrange P1 and P2 elements on triangles, scalar and vector
  valued, with Gaussian quadrature on intervals and triangles.
- Sparse assembly of stiffness, mass, advection and boundary operators
  and load vectors. Terms with coefficients that are constant on each
  element are assembled from reference-element matrices, exactly and
  without a quadrature loop.
- Coefficients given as numbers, functions of `(x, y)` or `(x, y, t)`,
  finite-element functions, constant or spatially varying symmetric
  tensors
  ([`tensor_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/tensor_coefficient.md))
  or piecewise values by region
  ([`region_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/region_coefficient.md)).
- Interpolation, exact transfer to refined meshes, and conservative L2
  projection between any two meshes of a domain
  ([`l2_project()`](https://adamrsophie.github.io/pdecompR/reference/l2_project.md)).
- Dirichlet, Neumann and Robin conditions; pure Neumann problems with a
  compatibility check and a zero-mean gauge on each connected component.
- Periodic conditions for elliptic problems
  ([`periodic_bc()`](https://adamrsophie.github.io/pdecompR/reference/periodic_bc.md)),
  by elimination of the identified degrees of freedom, including fully
  periodic problems with a zero-mean gauge.
- A weak-form language
  ([`trial_function()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
  [`test_function()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
  [`grad()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
  [`dot()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
  [`integrate_domain()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
  [`integrate_boundary()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md))
  for linear scalar problems.

### Problem classes

- Elliptic problems with diffusion, advection and reaction, with
  optional SUPG stabilization.
- Heat and general linear parabolic problems with backward Euler,
  Crank-Nicolson (Rannacher startup), BDF2 and L-stable SDIRK2, adaptive
  SDIRK2 stepping with an embedded error estimate, time-dependent data
  and coefficients, and mass diagnostics.
- Reaction-diffusion systems with IMEX and fully implicit Newton time
  stepping.
- Quasilinear diffusion, semilinear reaction and p-Laplace problems
  solved by damped Newton iteration.
- Linear elasticity in plane strain and plane stress with displacement,
  roller and traction conditions, and strain, stress and von Mises
  output.
- Damped wave equation and elastodynamics integrated by the Newmark
  average acceleration method, second order and energy conserving, with
  the discrete energy reported at every step.
- [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md)
  rejects elements that overlap without sharing an edge (folded node
  neighborhoods or crossing boundary edges).
- Polygon domains with holes (`pde_domain(holes = )`), meshed by
  bridging each hole to the outer boundary before ear clipping; hole
  sides carry their own markers.
- Incompressible Stokes flow with Taylor-Hood P2-P1 elements, velocity,
  pseudo-traction and do-nothing outflow conditions, solved by sparse LU
  or by GMRES with a block-diagonal preconditioner
  ([`stokes_problem()`](https://adamrsophie.github.io/pdecompR/reference/stokes_problem.md)).
- Eigenvalue problems for symmetric elliptic operators and elastic
  vibration modes by shift-and-invert subspace iteration
  ([`solve_eigen()`](https://adamrsophie.github.io/pdecompR/reference/solve_eigen.md)).

### Solvers

- Sparse Cholesky and LU with symmetric equilibration, and
  preconditioned conjugate gradients, BiCGSTAB and GMRES, with residual
  and conditioning diagnostics and an explicit status for every solve.
  Several right-hand sides are solved with one factorization, and time
  stepping warm starts iterative solvers.
- Jacobi, symmetric Gauss-Seidel (SSOR) and user-supplied
  preconditioners, and a geometric multigrid V-cycle on nested meshes
  ([`multigrid_preconditioner()`](https://adamrsophie.github.io/pdecompR/reference/multigrid_preconditioner.md))
  with mesh-independent iteration counts.

### Verification and adaptivity

- L2, H1 and energy error norms against exact or finite-element
  reference solutions on nested or unrelated meshes (integrated on the
  mesh overlay), and convergence studies with observed orders.
- Symbolic manufactured solutions through the optional ‘symcompR’
  package.
- A residual a posteriori error estimator, Doerfler marking and an
  adaptive solve loop.
- A benchmark suite of ten reference problems in `benchmarks/`.

### Inference and control

- Discrete forward and adjoint sensitivities of linear elliptic
  problems, with finite-difference or exact operator derivatives
  ([`affine_derivative()`](https://adamrsophie.github.io/pdecompR/reference/affine_derivative.md))
  and Taylor remainder checks
  ([`taylor_test()`](https://adamrsophie.github.io/pdecompR/reference/taylor_test.md)).
- Adjoint gradients of nonlinear problems (one transposed Newton
  Jacobian solve) and of backward Euler and Crank-Nicolson time stepping
  (a discrete backward sweep) through
  [`nonlinear_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/nonlinear_sensitivity.md)
  and
  [`transient_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/nonlinear_sensitivity.md).
- Parameter estimation with adjoint gradients, bounds (including fixed
  parameters), Gauss-Newton standard errors and identifiability
  diagnostics.
- Distributed optimal control with optional control bounds.
