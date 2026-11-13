# Boundary degrees of freedom

Boundary degrees of freedom

## Usage

``` r
boundary_dofs(space, marker = "boundary")
```

## Arguments

- space:

  An
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

- marker:

  Boundary markers; `"boundary"` selects the whole boundary.

## Value

Sorted integer indices of the degrees of freedom on the selected
boundary edges.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 2))
boundary_dofs(V, "bottom")
#> [1] 1 2 3
```
