# Subdomains (regions) of a mesh

Every triangle of a mesh belongs to a named region (subdomain); by
default all triangles belong to `"domain"`. `mesh_regions()` returns the
region of each triangle, `set_region()` assigns a region to the
triangles whose centroids satisfy a predicate, and
[`read_gmsh()`](https://adamrsophie.github.io/pdecompR/reference/read_gmsh.md)
takes regions from physical surfaces. Refinement keeps the region of the
parent triangle. Regions carry piecewise coefficients through
[`region_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/region_coefficient.md).

## Usage

``` r
mesh_regions(mesh)

set_region(mesh, region, which)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- region:

  Region name.

- which:

  A function `function(x, y)` returning `TRUE` for the centroids of the
  triangles to assign.

## Value

`mesh_regions()`: a character vector with one entry per triangle.
`set_region()`: the updated mesh.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 4)
mesh <- set_region(mesh, "inclusion", function(x, y) (x - 0.5)^2 + (y - 0.5)^2 < 0.1)
table(mesh_regions(mesh))
#> 
#>    domain inclusion 
#>        24         8 
```
