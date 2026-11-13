# Mesh geometry and quality

`mesh_quality()` returns per-element geometric measures: area, longest
and shortest edge, minimum and maximum interior angle (degrees), and the
normalized shape ratio `2 r / R` (inradius over circumradius, equal to 1
for an equilateral triangle and tending to 0 for degenerate ones).

## Usage

``` r
mesh_quality(mesh)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

## Value

A data frame with one row per element.

## Examples

``` r
q <- mesh_quality(pde_mesh_rectangle(nx = 2))
head(q)
#>    area     h_max h_min min_angle max_angle shape_ratio
#> 1 0.125 0.7071068   0.5        45        90   0.8284271
#> 2 0.125 0.7071068   0.5        45        90   0.8284271
#> 3 0.125 0.7071068   0.5        45        90   0.8284271
#> 4 0.125 0.7071068   0.5        45        90   0.8284271
#> 5 0.125 0.7071068   0.5        45        90   0.8284271
#> 6 0.125 0.7071068   0.5        45        90   0.8284271
```
