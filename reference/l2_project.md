# Conservative L2 projection

Computes the L2 projection `P f` of `f` onto `space`: the finite-element
function with `int (P f) v = int f v` for every `v` in the space.
Because constants belong to the space, the projection preserves the
integral of `f` over the whole domain, which interpolation does not. It
does not preserve integrals over individual elements or unions of them:
local conservation would need a discontinuous space or explicit
constraints.

## Usage

``` r
l2_project(f, space, quad_degree = NULL)
```

## Arguments

- f:

  Function to project.

- space:

  Target
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

- quad_degree:

  Quadrature exactness degree for the load integrals.

## Value

An
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
on `space`.

## Details

`f` may be a number, a `function(x, y)`, or an
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
on any mesh of the same domain. For a finite-element function on another
mesh the load integrals are computed exactly on the overlay of the two
meshes (pointwise on the target mesh when it refines the source mesh),
so the transfer is conservative up to rounding.

## Examples

``` r
a <- interpolate(function(x, y) sin(5 * x) * y, fe_space(pde_mesh_rectangle(nx = 12)))
V <- fe_space(pde_mesh_rectangle(nx = 5, pattern = "left"))
p <- l2_project(a, V)
c(sum(assemble_matrix(a$space, "mass") %*% a$values),
  sum(assemble_matrix(V, "mass") %*% p$values))
#> [1] 0.07003947 0.07003947
```
