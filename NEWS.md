# pdecompR 0.1.0

* First version. A pure R finite-element framework for two-dimensional
  scalar PDEs:
  * domains (rectangles, simple polygons), structured and polygon triangular
    meshes with validation, quality measures and named boundary markers,
    Gmsh 2.2 and 4.1 import, VTK export;
  * continuous Lagrange P1 and P2 elements, Gaussian quadrature on
    intervals and triangles, sparse assembly through 'Matrix';
  * Dirichlet, Neumann and Robin conditions, pure Neumann problems with
    compatibility checks and a zero-mean gauge;
  * linear elasticity with vector-valued P1 and P2 displacements, plane
    strain and plane stress, displacement (including roller) and traction
    conditions;
  * elliptic (diffusion, advection with optional SUPG, reaction), heat and
    general parabolic (backward Euler, Crank-Nicolson), reaction-diffusion
    systems (IMEX and Newton), quasilinear and p-Laplace solvers;
  * a weak-form language (`trial_function()`, `grad()`, `dot()`,
    `integrate_domain()`, `integrate_boundary()`);
  * error norms, convergence studies, a residual error estimator, Doerfler
    marking and newest vertex bisection;
  * discrete forward and adjoint sensitivities, parameter estimation and
    distributed optimal control.
