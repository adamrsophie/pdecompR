# Mesh refinement

`refine_uniform()` splits every triangle into four by joining its edge
midpoints (red refinement). Angles are preserved, so shape regularity is
inherited from the parent mesh.

## Usage

``` r
refine_uniform(mesh)

refine_mesh(mesh, marked)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- marked:

  Logical vector over elements, or integer element indices.

## Value

The refined `pde_mesh`.

## Details

`refine_mesh()` refines the marked elements by newest vertex bisection
and bisects neighbors as needed to keep the mesh conforming (no hanging
nodes). Every edge of a marked element is bisected, so marked elements
are split into four children. On the first call the longest edge of each
triangle is chosen as its refinement edge; afterwards the newest vertex
rule applies. Newest vertex bisection produces at most four similarity
classes of triangles per initial triangle, so the minimum angle stays
bounded away from zero.

Both functions preserve boundary markers: a bisected boundary edge
yields two edges with the parent's marker. The returned mesh records the
parent edge of each new node in `info$new_nodes` (a two-column matrix of
node indices), which makes exact transfer of P1 functions possible, see
[`transfer_function()`](https://adamrsophie.github.io/pdecompR/reference/transfer_function.md).

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 2)
nrow(refine_uniform(mesh)$triangles)
#> [1] 32
fine <- refine_mesh(mesh, marked = 1)
nrow(fine$triangles)
#> [1] 15
```
