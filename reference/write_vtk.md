# Write a VTK file

Writes a mesh, finite-element function or solution in the legacy ASCII
VTK unstructured-grid format, readable by ParaView and VisIt. Point data
are the vertex values; for P2 functions only the vertex values are
written (the triangles are linear cells).

## Usage

``` r
write_vtk(x, file, name = "u")
```

## Arguments

- x:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md),
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
  or `pde_solution`.

- file:

  Output path.

- name:

  Name of the point data array.

## Value

The file path, invisibly.

## Examples

``` r
f <- tempfile(fileext = ".vtk")
write_vtk(pde_mesh_rectangle(nx = 2), f)
```
