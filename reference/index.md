# Package index

## Domains and meshes

- [`pde_domain()`](https://adamrsophie.github.io/pdecompR/reference/pde_domain.md)
  : Computational domains
- [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md)
  : Triangular meshes
- [`pde_mesh_rectangle()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh_rectangle.md)
  : Structured triangular mesh of a rectangle
- [`pde_mesh_polygon()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh_polygon.md)
  : Triangular mesh of a polygon
- [`mesh_quality()`](https://adamrsophie.github.io/pdecompR/reference/mesh_quality.md)
  : Mesh geometry and quality
- [`mesh_size()`](https://adamrsophie.github.io/pdecompR/reference/mesh_size.md)
  : Mesh size
- [`boundary_markers()`](https://adamrsophie.github.io/pdecompR/reference/boundary_markers.md)
  [`boundary_nodes()`](https://adamrsophie.github.io/pdecompR/reference/boundary_markers.md)
  : Boundary queries
- [`set_boundary_marker()`](https://adamrsophie.github.io/pdecompR/reference/set_boundary_marker.md)
  : Relabel boundary edges
- [`mesh_regions()`](https://adamrsophie.github.io/pdecompR/reference/mesh_regions.md)
  [`set_region()`](https://adamrsophie.github.io/pdecompR/reference/mesh_regions.md)
  : Subdomains (regions) of a mesh
- [`region_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/region_coefficient.md)
  : Piecewise coefficients by region
- [`refine_uniform()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
  [`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
  : Mesh refinement
- [`read_gmsh()`](https://adamrsophie.github.io/pdecompR/reference/read_gmsh.md)
  : Read a Gmsh mesh
- [`write_vtk()`](https://adamrsophie.github.io/pdecompR/reference/write_vtk.md)
  : Write a VTK file

## Finite elements and quadrature

- [`fe_element()`](https://adamrsophie.github.io/pdecompR/reference/fe_element.md)
  : Lagrange finite elements on triangles
- [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md)
  : Finite-element function spaces
- [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
  [`interpolate()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
  : Finite-element functions
- [`boundary_dofs()`](https://adamrsophie.github.io/pdecompR/reference/boundary_dofs.md)
  : Boundary degrees of freedom
- [`evaluate()`](https://adamrsophie.github.io/pdecompR/reference/evaluate.md)
  [`evaluation_matrix()`](https://adamrsophie.github.io/pdecompR/reference/evaluate.md)
  : Evaluate finite-element functions at points
- [`transfer_function()`](https://adamrsophie.github.io/pdecompR/reference/transfer_function.md)
  : Transfer a finite-element function to another mesh
- [`quadrature_interval()`](https://adamrsophie.github.io/pdecompR/reference/quadrature_interval.md)
  [`quadrature_triangle()`](https://adamrsophie.github.io/pdecompR/reference/quadrature_interval.md)
  : Numerical quadrature rules
- [`integrate_mesh()`](https://adamrsophie.github.io/pdecompR/reference/integrate_mesh.md)
  : Integrate a function over a mesh
- [`element_stiffness_p1()`](https://adamrsophie.github.io/pdecompR/reference/element_stiffness_p1.md)
  [`element_mass_p1()`](https://adamrsophie.github.io/pdecompR/reference/element_stiffness_p1.md)
  : Local P1 element matrices

## Assembly and linear algebra

- [`assemble_matrix()`](https://adamrsophie.github.io/pdecompR/reference/assemble_matrix.md)
  [`assemble_vector()`](https://adamrsophie.github.io/pdecompR/reference/assemble_matrix.md)
  [`assemble_boundary_matrix()`](https://adamrsophie.github.io/pdecompR/reference/assemble_matrix.md)
  [`assemble_boundary_vector()`](https://adamrsophie.github.io/pdecompR/reference/assemble_matrix.md)
  : Finite-element assembly
- [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
  [`linear_solver()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
  : Sparse linear solvers
- [`multigrid_preconditioner()`](https://adamrsophie.github.io/pdecompR/reference/multigrid_preconditioner.md)
  : Geometric multigrid preconditioner

## Problems and solvers

- [`dirichlet_bc()`](https://adamrsophie.github.io/pdecompR/reference/dirichlet_bc.md)
  [`neumann_bc()`](https://adamrsophie.github.io/pdecompR/reference/dirichlet_bc.md)
  [`robin_bc()`](https://adamrsophie.github.io/pdecompR/reference/dirichlet_bc.md)
  : Boundary conditions
- [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md)
  [`poisson_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md)
  : Linear elliptic problems
- [`parabolic_problem()`](https://adamrsophie.github.io/pdecompR/reference/parabolic_problem.md)
  [`heat_problem()`](https://adamrsophie.github.io/pdecompR/reference/parabolic_problem.md)
  : Linear parabolic problems
- [`reaction_diffusion_problem()`](https://adamrsophie.github.io/pdecompR/reference/reaction_diffusion_problem.md)
  : Reaction-diffusion systems
- [`nonlinear_problem()`](https://adamrsophie.github.io/pdecompR/reference/nonlinear_problem.md)
  [`p_laplace_problem()`](https://adamrsophie.github.io/pdecompR/reference/nonlinear_problem.md)
  : Nonlinear stationary problems
- [`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md)
  : Linear elasticity
- [`elasticity_stress()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_stress.md)
  : Stresses and strains of an elasticity solution
- [`displacement_bc()`](https://adamrsophie.github.io/pdecompR/reference/displacement_bc.md)
  [`traction_bc()`](https://adamrsophie.github.io/pdecompR/reference/displacement_bc.md)
  : Boundary conditions for linear elasticity
- [`trial_function()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  [`test_function()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  [`grad()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  [`vector_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  [`dot()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  [`integrate_domain()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  [`integrate_boundary()`](https://adamrsophie.github.io/pdecompR/reference/trial_function.md)
  : Weak-form language
- [`variational_problem()`](https://adamrsophie.github.io/pdecompR/reference/variational_problem.md)
  : Variational problems
- [`discretize()`](https://adamrsophie.github.io/pdecompR/reference/discretize.md)
  : Discretize a PDE problem
- [`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
  : Solve a PDE problem
- [`update_mesh()`](https://adamrsophie.github.io/pdecompR/reference/update_mesh.md)
  : Replace the mesh of a problem
- [`solution_at()`](https://adamrsophie.github.io/pdecompR/reference/solution_at.md)
  : Extract a time snapshot

## Errors and adaptivity

- [`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)
  : Errors against a reference solution
- [`convergence_study()`](https://adamrsophie.github.io/pdecompR/reference/convergence_study.md)
  : Mesh convergence studies
- [`manufactured_solution()`](https://adamrsophie.github.io/pdecompR/reference/manufactured_solution.md)
  : Manufactured solutions with symbolic derivatives
- [`estimate_error()`](https://adamrsophie.github.io/pdecompR/reference/estimate_error.md)
  : A posteriori error estimation
- [`mark_elements()`](https://adamrsophie.github.io/pdecompR/reference/mark_elements.md)
  : Mark elements for refinement
- [`adaptive_solve()`](https://adamrsophie.github.io/pdecompR/reference/adaptive_solve.md)
  : Adaptive finite-element solution

## Sensitivities, inverse problems and control

- [`pde_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/pde_sensitivity.md)
  : Parameter sensitivities of linear elliptic problems
- [`affine_derivative()`](https://adamrsophie.github.io/pdecompR/reference/affine_derivative.md)
  : Exact derivatives of an affine problem family
- [`taylor_test()`](https://adamrsophie.github.io/pdecompR/reference/taylor_test.md)
  : Taylor remainder test of an objective gradient
- [`pde_observations()`](https://adamrsophie.github.io/pdecompR/reference/pde_observations.md)
  [`misfit_objective()`](https://adamrsophie.github.io/pdecompR/reference/pde_observations.md)
  : Point observations and data misfit
- [`estimate_parameters()`](https://adamrsophie.github.io/pdecompR/reference/estimate_parameters.md)
  : Parameter estimation for PDE models
- [`optimal_control()`](https://adamrsophie.github.io/pdecompR/reference/optimal_control.md)
  : Distributed optimal control

## Visualization

- [`plot(`*`<pde_vector_solution>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot(`*`<pde_mesh>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot(`*`<fe_function>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot(`*`<pde_solution>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot(`*`<pde_transient>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot(`*`<pde_convergence>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot_indicators()`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  [`plot(`*`<pde_adaptive>`*`)`](https://adamrsophie.github.io/pdecompR/reference/pde_plots.md)
  : Plot meshes and finite-element functions

## Package

- [`pdecompR`](https://adamrsophie.github.io/pdecompR/reference/pdecompR-package.md)
  [`pdecompR-package`](https://adamrsophie.github.io/pdecompR/reference/pdecompR-package.md)
  : pdecompR: Partial Differential Equations and Finite-Element
  Computing in R
