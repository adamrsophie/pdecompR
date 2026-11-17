# Computational domains

`pde_domain()` describes the geometry on which a PDE is posed. Domains
are pure geometric descriptions; a mesh is built from a domain with
[`pde_mesh_rectangle()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh_rectangle.md)
or
[`pde_mesh_polygon()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh_polygon.md).
One-dimensional and three-dimensional domains are planned but not yet
available.

## Usage

``` r
pde_domain(
  type = c("rectangle", "polygon"),
  x = c(0, 1),
  y = c(0, 1),
  vertices = NULL,
  markers = "boundary",
  holes = NULL,
  hole_markers = "hole"
)
```

## Arguments

- type:

  One of `"rectangle"` or `"polygon"`.

- x:

  For rectangles, the two endpoints in the x direction.

- y:

  For rectangles, the two endpoints in the y direction.

- vertices:

  For polygons, a two-column numeric matrix of vertices.

- markers:

  For polygons, a character vector naming each side, or a single name
  used for every side.

- holes:

  For polygons, an optional list of vertex matrices, each a simple
  polygon strictly inside the outer polygon; holes must not touch each
  other or the outer boundary.

- hole_markers:

  Marker of the sides of each hole: one name for all holes or one per
  hole.

## Value

An object of class `pde_domain`.

## Details

Supported types:

- `"rectangle"`: the open rectangle `(x[1], x[2]) x (y[1], y[2])`.
  Boundary markers are `"bottom"`, `"right"`, `"top"` and `"left"`.

- `"polygon"`: a simple polygon given by its vertices in order. Side `i`
  joins vertex `i` to vertex `i + 1` (and the last vertex to the first).
  `markers` names each side; by default all sides are `"boundary"`.

Every domain also answers to the marker `"boundary"`, which selects the
whole boundary.

Validation rejects non-finite coordinates, zero-length sides, repeated
vertices, zero area and self-intersecting polygons. A clockwise polygon
is reversed to counterclockwise order (its side markers are reversed
with it) and a message reports the change.

## Examples

``` r
pde_domain("rectangle", x = c(0, 1), y = c(0, 2))
#> <pde_domain> rectangle 
#>   (0, 1) x (0, 2)
#>   markers: bottom, right, top, left 
L <- rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
pde_domain("polygon", vertices = L)
#> <pde_domain> polygon 
#>   vertices: 6 
#>   markers: boundary 
```
