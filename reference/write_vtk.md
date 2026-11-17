# Write a VTK file

Writes a mesh, finite-element function or solution in the legacy ASCII
VTK unstructured-grid format, readable by ParaView and VisIt. P1 data
are written on linear triangles; P2 data on quadratic triangles (VTK
cell type 22) with all six nodal values, which ParaView renders with
curved interpolation. A displacement solution of
[`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md)
is written as a `VECTORS` array, optionally with the nodal stresses of
[`elasticity_stress()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_stress.md)
(on the vertices, as P1 data). A time-dependent solution is written at
one stored time (all species for reaction-diffusion systems); call it
once per time to write a series.

## Usage

``` r
write_vtk(x, file, name = "u", time = NULL, stress = FALSE)
```

## Arguments

- x:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md),
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md),
  `pde_solution`, `pde_vector_solution` or `pde_transient`.

- file:

  Output path.

- name:

  Name of the point data array (scalar data and displacements).

- time:

  Stored time to write for a `pde_transient` (default: the last).

- stress:

  For elasticity solutions, also write `sxx`, `syy`, `sxy` and
  `von_mises` (P1 solutions only, since the stresses are vertex values).

## Value

The file path, invisibly.

## Examples

``` r
f <- tempfile(fileext = ".vtk")
write_vtk(pde_mesh_rectangle(nx = 2), f)
```
