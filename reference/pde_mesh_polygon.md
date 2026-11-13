# Triangular mesh of a polygon

Triangulates a simple polygon by ear clipping, improves the coarse
triangulation with Lawson edge flips (making it Delaunay while keeping
the polygon sides), and then refines it uniformly until the longest edge
is at most `h`. Uniform refinement keeps the angles of the coarse
triangulation, so mesh quality is set by the polygon's shape.

## Usage

``` r
pde_mesh_polygon(domain, h = Inf, markers = "boundary")
```

## Arguments

- domain:

  A polygon or rectangle
  [`pde_domain()`](https://adamrsophie.github.io/pdecompR/reference/pde_domain.md),
  or a vertex matrix.

- h:

  Target maximum edge length.

- markers:

  Side markers when `domain` is a vertex matrix.

## Value

A `pde_mesh` whose boundary edges carry the polygon side markers.

## Details

This is a simple, dependency-free mesher for polygonal domains. It is
not a quality-guaranteeing constrained Delaunay mesher; for complex
geometry, generate the mesh with Gmsh and read it with
[`read_gmsh()`](https://adamrsophie.github.io/pdecompR/reference/read_gmsh.md).

## Examples

``` r
L <- rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
mesh <- pde_mesh_polygon(L, h = 0.5)
summary(mesh)
#> Triangular mesh summary
#>   nodes 45  elements 64  edges 108  boundary edges 24 
#>   total area         3 
#>   edge length: min  0.25  max   0.5 
#>   angles (deg): min    45  max    90 
#>   shape ratio 2r/R: min 0.8284  median 0.8284  max 0.8284 
#>   boundary edges by marker:
#>     boundary: 24
```
