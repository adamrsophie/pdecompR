# Boundary queries

`boundary_markers()` lists the markers present on a mesh.
`boundary_nodes()` returns the sorted indices of nodes on boundary edges
with the given markers. Selection is by marker only; no floating-point
coordinate comparisons are involved.

## Usage

``` r
boundary_markers(mesh)

boundary_nodes(mesh, marker = "boundary")
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- marker:

  Character vector of markers. `"boundary"` selects every boundary edge.

## Value

A character vector of markers, or an integer vector of node indices.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 2)
boundary_markers(mesh)
#> [1] "bottom" "right"  "top"    "left"  
boundary_nodes(mesh, "left")
#> [1] 1 4 7
```
