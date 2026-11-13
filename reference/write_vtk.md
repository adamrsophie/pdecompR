# Write a VTK file

Writes a mesh, finite-element function or solution in the legacy ASCII
VTK unstructured-grid format, readable by ParaView and VisIt. Point data
are the vertex values; for P2 functions only the vertex values are
written (the triangles are linear cells). A displacement solution of
[`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md)
is written as a `VECTORS` array; a time-dependent solution is written at
one stored time (all species for reaction-diffusion systems). Call it
once per time to write a series.

## Usage

``` r
write_vtk(x, file, name = "u", time = NULL)
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

## Value

The file path, invisibly.

## Examples

``` r
f <- tempfile(fileext = ".vtk")
write_vtk(pde_mesh_rectangle(nx = 2), f)
```
