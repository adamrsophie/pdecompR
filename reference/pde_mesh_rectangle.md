# Structured triangular mesh of a rectangle

Splits the rectangle into `nx` by `ny` cells and each cell into two
triangles.

## Usage

``` r
pde_mesh_rectangle(
  x = c(0, 1),
  y = c(0, 1),
  nx = 8,
  ny = nx,
  pattern = c("right", "left", "alternating")
)
```

## Arguments

- x, y:

  Endpoints of the rectangle in each direction.

- nx, ny:

  Number of cells in each direction.

- pattern:

  Diagonal pattern: `"right"` (diagonals from lower left to upper
  right), `"left"` (lower right to upper left) or `"alternating"` (the
  diagonal direction alternates from cell to cell).

## Value

A `pde_mesh` with boundary markers `"bottom"`, `"right"`, `"top"` and
`"left"`. Corner nodes belong to both adjacent sides.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 4, ny = 4)
summary(mesh)
#> Triangular mesh summary
#>   nodes 25  elements 32  edges 56  boundary edges 16 
#>   total area         1 
#>   edge length: min  0.25  max 0.3536 
#>   angles (deg): min    45  max    90 
#>   shape ratio 2r/R: min 0.8284  median 0.8284  max 0.8284 
#>   boundary edges by marker:
#>     bottom: 4
#>     left: 4
#>     right: 4
#>     top: 4
```
