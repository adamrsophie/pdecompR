# Relabel boundary edges

Assigns a new marker to the boundary edges whose midpoints satisfy a
predicate. This is a convenience for imported meshes without physical
markers; selection uses edge midpoints, so it is robust for straight
sides as long as `which` is not evaluated exactly on a corner.

## Usage

``` r
set_boundary_marker(mesh, marker, which)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- marker:

  New marker name.

- which:

  A function `function(x, y)` returning `TRUE` for edge midpoints to
  relabel.

## Value

The mesh with updated markers.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 4)
mesh <- set_boundary_marker(mesh, "inlet", function(x, y) x < 1e-8 & y > 0.5)
table(mesh$boundary$marker)
#> 
#> bottom  inlet   left  right    top 
#>      4      2      2      4      4 
```
