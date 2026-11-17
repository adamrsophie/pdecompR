# Curved boundaries

Declares the boundary edges with marker `marker` as approximations of a
curve, given by its closest point `projection(x, y)`. Later calls of
[`refine_uniform()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
and
[`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
then place new boundary nodes on the curve instead of at chord
midpoints, so the polygonal domain converges to the curved one and P1
elements keep their second-order L2 accuracy (a mesh that keeps the
coarse chords stalls at the geometric error). The elements stay
straight-sided: P2 edge nodes lie on the chords, which limits P2 to the
accuracy of the polygonal domain (isoparametric elements are not
provided).

## Usage

``` r
curved_boundary(mesh, marker, projection)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- marker:

  Boundary marker of the curved edges.

- projection:

  `function(x, y)` returning a two-column matrix of the closest points
  on the curve.

## Value

The mesh, with the curve stored in `curves`.

## Details

The existing nodes on the marked edges must already lie on the curve.
The projection is kept by refinement and by meshes derived from the
refined ones; relabeling edges with
[`set_boundary_marker()`](https://adamrsophie.github.io/pdecompR/reference/set_boundary_marker.md)
to a new marker drops the curve for them. Because refined meshes no
longer cover the same domain as their parents,
[`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)
against a finite-element reference from another mesh of the sequence and
[`l2_project()`](https://adamrsophie.github.io/pdecompR/reference/l2_project.md)
between them stop with an error; compare with exact functions instead.
Multigrid prolongation and
[`transfer_function()`](https://adamrsophie.github.io/pdecompR/reference/transfer_function.md)
use the chord average at projected nodes.

## Examples

``` r
n <- 16
circle <- cbind(cos(2 * pi * (1:n) / n), sin(2 * pi * (1:n) / n))
disk <- curved_boundary(pde_mesh_polygon(circle), "boundary",
                        function(x, y) cbind(x, y) / sqrt(x^2 + y^2))
fine <- refine_uniform(refine_uniform(disk))
range(sqrt(rowSums(fine$nodes[boundary_nodes(fine), ]^2)))
#> [1] 1 1
```
