# Finite-element function spaces

`fe_space()` builds a continuous Lagrange space on a triangular mesh.

## Usage

``` r
fe_space(mesh, degree = 1, family = "Lagrange")
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- degree:

  Polynomial degree, 1 or 2.

- family:

  Element family. Only `"Lagrange"` is implemented.

## Value

An object of class `fe_space` with fields `mesh`, `element`, `ndof`,
`cell_dofs` and `dof_coords`.

## Degrees of freedom

For P1 the degrees of freedom are the mesh nodes, in node order. For P2
they are the nodes followed by the edge midpoints: the midpoint of mesh
edge `k` has index `nrow(mesh$nodes) + k`. `cell_dofs` maps each element
to its global degrees of freedom in local basis order.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 4), degree = 2)
V
#> <fe_space> continuous Lagrange P2, 81 degrees of freedom on 32 triangles
```
