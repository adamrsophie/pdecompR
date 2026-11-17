# Triangular meshes

`pde_mesh()` builds and validates a two-dimensional triangular mesh from
node coordinates and element connectivity.

## Usage

``` r
pde_mesh(
  nodes,
  triangles,
  boundary_edges = NULL,
  boundary_markers = NULL,
  fix_orientation = FALSE,
  regions = NULL,
  check_overlap = TRUE
)
```

## Arguments

- nodes:

  Numeric matrix with two columns (x, y).

- triangles:

  Integer matrix with three columns of node indices.

- boundary_edges:

  Optional two-column matrix of node pairs on the boundary. When given,
  it must list every boundary edge exactly once.

- boundary_markers:

  Character vector with one marker per row of `boundary_edges`. Ignored
  when `boundary_edges` is `NULL`; all boundary edges are then marked
  `"boundary"`.

- fix_orientation:

  Reorder clockwise triangles instead of failing.

- regions:

  Optional character vector with one region (subdomain) name per
  triangle, for piecewise coefficients
  ([`region_coefficient()`](https://adamrsophie.github.io/pdecompR/reference/region_coefficient.md)).
  Default: every triangle belongs to region `"domain"`.

- check_overlap:

  Check that elements do not overlap geometrically (see the Validation
  section). Turn it off for meshes that intentionally contain coincident
  nodes, such as slits or cracks with duplicated nodes.

## Value

An object of class `pde_mesh` with components `nodes`, `triangles`,
`edges`, `element_edges`, `edge_elements`, `boundary` (a data frame of
boundary edges) and `info`.

## Conventions

- Indices are one-based, as everywhere in R. Node `i` is row `i` of
  `nodes`; element `e` is row `e` of `triangles`.

- Triangles are stored counterclockwise (positive signed area).

- The local edges of a triangle `(v1, v2, v3)` are edge 1 `(v1, v2)`,
  edge 2 `(v2, v3)` and edge 3 `(v3, v1)`.

- Boundary edges are oriented so that the domain lies on their left;
  their outward normal is the edge direction rotated clockwise.

- Each boundary edge carries exactly one marker (a character label). The
  label `"boundary"` always selects all boundary edges.

## Validation

The constructor rejects non-finite coordinates, invalid or repeated node
indices, duplicate elements, degenerate (zero-area) triangles, unused
nodes, edges shared by more than two triangles, and boundary definitions
that do not match the topological boundary. Clockwise triangles are
rejected unless `fix_orientation = TRUE`, in which case they are
reordered and the number of reordered elements is reported.

With `check_overlap = TRUE` it also rejects elements that overlap
without sharing an edge: the corner angles around every interior node
must sum to `2 pi` (at most `2 pi` at boundary nodes) and boundary edges
without a common node must not touch.

## Examples

``` r
nodes <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
tris <- rbind(c(1, 2, 3), c(1, 3, 4))
mesh <- pde_mesh(nodes, tris)
mesh
#> <pde_mesh> 2D triangular mesh
#>   nodes: 4  elements: 2  edges: 5  boundary edges: 4 
#>   h = 1.414 
#>   markers: boundary 
```
