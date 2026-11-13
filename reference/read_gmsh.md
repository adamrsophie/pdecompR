# Read a Gmsh mesh

Reads two-dimensional triangular meshes in the ASCII Gmsh formats 2.2
and 4.1. Triangles (element type 2) become mesh elements; line elements
(type 1) on the boundary carry boundary markers taken from their
physical group names (or numbers when names are absent). Other element
types (points, quadrilaterals, higher-order elements) are ignored with a
message. Boundary edges not covered by a physical line receive the
marker `"boundary"`. Nodes not used by triangles are dropped and the
remaining nodes renumbered; clockwise triangles are reoriented.

## Usage

``` r
read_gmsh(file)
```

## Arguments

- file:

  Path to a `.msh` file.

## Value

A
[`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

## Examples

``` r
f <- system.file("extdata", "square.msh", package = "pdecompR")
mesh <- read_gmsh(f)
#> Ignored 1 element(s) that are not triangles or lines.
table(mesh$boundary$marker)
#> 
#> bottom   left  right    top 
#>      2      2      2      2 
```
