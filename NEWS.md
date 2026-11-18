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
  quadratic P2 cells, displacement vectors and stress fields, in legacy and
  XML (`.vtu`, raw binary or text) formats, and ParaView time series
  (`write_pvd()`); field and file names are XML escaped.

## Discretization

* Continuous Lagrange P1 and P2 elements on triangles, scalar and vector
  valued, with Gaussian quadrature on intervals and triangles.
* Sparse assembly of stiffness, mass, advection and boundary operators and
  load vectors. Terms with coefficients that are constant on each element
  are assembled from reference-element matrices, exactly and without a
  quadrature loop.
* Coefficients given as numbers, functions of `(x, y)` or `(x, y, t)`,
  finite-element functions, constant or spatially varying symmetric tensors
  (`tensor_coefficient()`) or piecewise values by region
  (`region_coefficient()`).
* Interpolation, exact transfer to refined meshes, and L2 projection
  between any two meshes of the same domain (`l2_project()`), which
  preserves the integral over the whole domain. Meshes of different
  domains (for example a subdomain, or a refined curved boundary) are
  rejected.
* Dirichlet, Neumann and Robin conditions; pure Neumann problems with a
  compatibility check and a zero-mean gauge on each connected component.
  Data whose discrete compatibility defect exceeds 1 percent are rejected
  on every mesh; the gauge reports whether the original problem was solved.
* Periodic conditions for elliptic problems (`periodic_bc()`), by
  elimination of the identified degrees of freedom, including fully
  periodic problems with a zero-mean gauge. Matching is one to one in
  coordinates relative to the boundary spacing, so it does not depend on
  units or origin; prescribed values must agree within each periodic class.
  The same constraint is applied by `solve_eigen()`, `pde_sensitivity()`,
  `estimate_parameters()` and `optimal_control()`.
* A weak-form language (`trial_function()`, `test_function()`, `grad()`,
  `dot()`, `integrate_domain()`, `integrate_boundary()`) for linear scalar
  problems.

## Problem classes

* Elliptic problems with diffusion, advection and reaction, with optional
  SUPG stabilization.
* Heat and general linear parabolic problems with backward Euler,
  Crank-Nicolson (Rannacher startup), BDF2 and L-stable SDIRK2, adaptive
  SDIRK2 stepping with an embedded error estimate, time-dependent data and
  coefficients, and mass diagnostics.
* Reaction-diffusion systems with IMEX and fully implicit Newton time
  stepping.
* Quasilinear diffusion, semilinear reaction and p-Laplace problems solved
  by damped Newton iteration. Without Dirichlet data, a converged state
  that stays a solution when shifted by a constant (large or small shifts,
  so a flat interval of solutions is found too) is reported as not unique.
  Newton convergence floors are relative to the size of the terms of the
  equations, so they do not depend on units, and the default initial guess
  balances a reaction that vanishes to first order at zero.
* Linear elasticity in plane strain and plane stress with displacement,
  roller and traction conditions, and strain, stress and von Mises output.
* Damped wave equation and elastodynamics integrated by the Newmark
  average acceleration method, second order and energy conserving, with
  the discrete energy reported at every step.
* `pde_mesh()` rejects elements that overlap without sharing an edge
  (folded node neighborhoods, crossing boundary edges, or a disconnected
  part lying inside another component, including on an interior edge).
  Refinement that moves nodes onto curves validates the result.
* `curved_boundary()` declares boundary parts as curves; refinement then
  places new boundary nodes on the curve, which restores second-order P1
  convergence on curved domains. Elements stay straight sided, so P2 is
  limited by the polygonal geometry (H1 rate 1.5).
* `amg_preconditioner()` adds smoothed aggregation algebraic multigrid for
  CG, with nearly constant iteration counts on unstructured and P2 meshes
  where geometric multigrid does not apply.
* `elasticity_problem(formulation = "mixed")` solves nearly incompressible
  elasticity with P2 displacements and a P1 pressure, free of volumetric
  locking; `elasticity_stress()` uses the pressure unknown.
* Every time stepper accepts a per-step `callback` and can `restart` from
  the final state of an earlier run (including multistep history, Newmark
  velocity and acceleration, the adaptive step size and startup steps), so
  long runs can be streamed and checkpointed. The state records its
  finite-element space (degree, dof coordinates and element connectivity),
  and restarts on another mesh, degree or numbering are rejected;
  reaction-diffusion species are matched by name and the initial condition
  is not evaluated on restart. Newmark restarts recompute the free
  acceleration from the loads of the new run.
* Time grids are built in elapsed time and end exactly at the requested
  time; steps that cannot be told apart at the time offset are rejected
  instead of producing repeated time stamps. Time tolerances are a few
  units in the last place, so large time origins work. Newmark startup
  never evaluates boundary data after the first step. Tensor coefficients
  with a time-dependent entry are reassembled in time, and adaptive SDIRK2
  under SUPG measures errors in the Galerkin mass norm.
* Time-dependent solvers report the range of nodal values over all steps
  (`diagnostics$range`) next to the integral, to expose violations of
  positivity or of the discrete maximum principle.
* Polygon domains with holes (`pde_domain(holes = )`), meshed by bridging
  each hole to the outer boundary before ear clipping; hole sides carry
  their own markers.
* Incompressible Stokes flow with Taylor-Hood P2-P1 elements, velocity,
  pseudo-traction and do-nothing outflow conditions, solved by sparse LU or
  by GMRES with a block-diagonal preconditioner (`stokes_problem()`). The
  pressure gauge is found from the constrained divergence operator per
  pressure component, each against its own coupling scale, so walls with
  only the normal velocity prescribed are gauged too. Free velocity
  translations and velocity spaces too small to control the pressure are
  reported.
* Eigenvalue problems for symmetric elliptic operators and elastic
  vibration modes by shift-and-invert subspace iteration
  (`solve_eigen()`).

## Solvers

* Sparse Cholesky and LU with symmetric equilibration, and preconditioned
  conjugate gradients, BiCGSTAB and GMRES, with residual and conditioning
  diagnostics and an explicit status for every solve. Several right-hand
  sides are solved with one factorization, and time stepping warm starts
  iterative solvers.
* Jacobi, symmetric Gauss-Seidel (SSOR) and user-supplied preconditioners,
  and a geometric multigrid V-cycle on nested meshes
  (`multigrid_preconditioner()`) with mesh-independent iteration counts;
  coarse levels inherit the homogeneous Dirichlet constraints. Algebraic
  multigrid is invariant to scaling of the system. GMRES is flexible, so
  preconditioners that are not fixed linear maps are applied consistently.
* The symmetry test behind the automatic Cholesky choice works on the
  equilibrated matrix, so it does not depend on the scale of the system.
* Pure Neumann gauges are decided per connected component (joined by
  periodic conditions) after Dirichlet, reaction and Robin terms; anchored
  components are solved as usual. The compatibility test normalizes the
  load and `compatibility_tol` is at most 0.01.
* Repeated assembly (time-dependent coefficients, Newton iterations)
  builds its sparse plans only once enough fills amortize them, and load
  vectors are not reassembled when only Dirichlet data move in time.

## Verification and adaptivity

* L2, H1 and energy error norms against exact or finite-element reference
  solutions on nested or unrelated meshes of the same domain (integrated on
  the mesh overlay, which keeps arbitrarily small cells), and convergence
  studies with observed orders. Common domains are certified triangle by
  triangle, so a missing small component is rejected, and the overlay works
  in a translated local frame. A reference on the same space is compared
  through the coefficient difference. The finite-difference gradient used
  when no exact gradient is given does not change when the domain is
  translated.
* Symbolic manufactured solutions through the optional 'symcompR' package.
* A residual a posteriori error estimator, Doerfler marking and an adaptive
  solve loop.
* A benchmark suite of ten reference problems in `benchmarks/`.

## Inference and control

* Discrete forward and adjoint sensitivities of linear elliptic problems,
  with finite-difference or exact operator derivatives
  (`affine_derivative()`) and Taylor remainder checks (`taylor_test()`).
* Adjoint gradients of nonlinear problems (one transposed Newton Jacobian
  solve) and of backward Euler and Crank-Nicolson time stepping (a
  discrete backward sweep) through `nonlinear_sensitivity()` and
  `transient_sensitivity()`; the transient sweep accumulates the gradient
  step by step, includes the derivative of the SUPG mass matrix, and can
  checkpoint the forward states. It keeps its own record of every step, so
  `output_times` only thins the returned solution, and it checks that the
  Dirichlet data do not depend on the parameters at every time level.
* Parameter estimation with adjoint gradients, bounds (including fixed
  parameters), Gauss-Newton standard errors and identifiability
  diagnostics: the rank of the observation sensitivities against the
  number of free parameters, and the null space of unidentified
  directions. The regularization of fixed parameters stays in the reported
  objective.
* Distributed optimal control with optional control bounds; contradictory
  bounds are rejected.
