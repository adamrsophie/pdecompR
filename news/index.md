# Changelog

## pdecompR 0.1.0

- First version. A pure R finite-element framework for two-dimensional
  scalar PDEs:
  - domains (rectangles, simple polygons), structured and polygon
    triangular meshes with validation, quality measures and named
    boundary markers, Gmsh 2.2 and 4.1 import, VTK export;
  - continuous Lagrange P1 and P2 elements, Gaussian quadrature on
    intervals and triangles, sparse assembly through ‘Matrix’;
  - Dirichlet, Neumann and Robin conditions, pure Neumann problems with
    compatibility checks and a zero-mean gauge;
  - linear elasticity with vector-valued P1 and P2 displacements, plane
    strain and plane stress, displacement (including roller) and
    traction conditions;
  - elliptic (diffusion, advection with optional SUPG, reaction), heat
    and general parabolic (backward Euler, Crank-Nicolson),
    reaction-diffusion systems (IMEX and Newton), quasilinear and
    p-Laplace solvers;
  - a weak-form language
    ([`trial_function()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
    [`grad()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
    [`dot()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
    [`integrate_domain()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md),
    [`integrate_boundary()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md));
  - error norms, convergence studies, a residual error estimator,
    Doerfler marking and newest vertex bisection;
  - discrete forward and adjoint sensitivities, parameter estimation and
    distributed optimal control.
- Robustness before release:
  - singular systems are detected from the assembled operator: Robin
    conditions with `alpha = 0`, variational problems, advection with
    zero-flux boundaries (left null vector compatibility) and elasticity
    constraints that leave a rigid motion;
  - linear solves that do not reach their tolerance give a warning, and
    iterative solvers guard against breakdown and zero right-hand sides;
  - the error estimator evaluates the diffusion on each side of an edge,
    gives Dirichlet conditions precedence as the solver does, and
    measures how well the space represents Dirichlet data;
  - finite-element coefficients and initial data are evaluated on their
    own space, so they may come from another mesh or degree;
    [`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)
    accepts a reference `fe_function`;
  - time-dependent Neumann and Robin data work in heat and
    reaction-diffusion problems; only a third argument named `t` makes a
    function time dependent;
  - Crank-Nicolson starts with two backward Euler steps (`startup`) to
    damp oscillations from rough initial data;
  - the weak-form language rejects products of directional derivatives,
    gradients of coefficient products and mixed function spaces;
  - SUPG accepts tensor diffusion; P2 functions are plotted on
    subdivided elements; plots compose with `par(mfrow = ...)`;
    [`write_vtk()`](https://adamrsophie.github.io/pdecompR/reference/write_vtk.md)
    writes displacement vectors and time-dependent solutions;
  - stricter input checks for reaction-diffusion species, boundary
    edges, observations, Dirichlet values, polygon size and Gmsh files.

### Additional features

- [`manufactured_solution()`](https://adamrsophie.github.io/pdecompR/reference/manufactured_solution.md)
  derives source terms and gradients symbolically with ‘symcompR’
  (optional) for verification studies.
- Mesh regions: `pde_mesh(regions = )`,
  [`set_region()`](https://adamrsophie.github.io/pdecompR/reference/mesh_regions.md),
  [`mesh_regions()`](https://adamrsophie.github.io/pdecompR/reference/mesh_regions.md),
  Gmsh physical surfaces, and
  [`region_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/region_coefficient.md)
  for piecewise coefficients that keep their geometry under refinement.
- [`elasticity_stress()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_stress.md)
  returns strains, stresses and the von Mises stress per element or
  averaged to vertices; `write_vtk(stress = TRUE)` writes them.
- [`write_vtk()`](https://adamrsophie.github.io/pdecompR/reference/write_vtk.md)
  writes P2 data on quadratic triangles (VTK cell type 22).
- Parabolic problems accept time-dependent diffusion, advection and
  reaction coefficients (reassembled each step, second order with
  Crank-Nicolson).
- Heat and reaction-diffusion solutions record the integral of each
  species at the stored times (`diagnostics$integral`), shown by
  [`summary()`](https://rdrr.io/r/base/summary.html).
- Direct solvers report `"singular"` when the factorization has a pivot
  ratio near machine precision; nonlinear problems without Dirichlet,
  Robin or reaction terms stop with a message.
- `benchmarks/reference-problems.R` runs the ten reference problems and
  records errors, observed orders, timings and peak memory.
