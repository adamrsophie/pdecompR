# Evaluate finite-element functions at points

`evaluate()` returns values (and optionally gradients) of a
finite-element function at arbitrary points. `evaluation_matrix()`
returns the sparse matrix `P` with `P %*% values` equal to the point
values, which is the observation operator used for inverse problems.

## Usage

``` r
evaluate(u, x, y = NULL, gradient = FALSE)

evaluation_matrix(space, x, y = NULL)
```

## Arguments

- u:

  An
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md).

- x, y:

  Coordinates (vectors of equal length), or a two-column matrix in `x`.

- gradient:

  Also return the gradient.

- space:

  An
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

## Value

`evaluate()` returns a numeric vector, or with `gradient = TRUE` a data
frame with columns `value`, `dx`, `dy`. `evaluation_matrix()` returns a
sparse `dgCMatrix`.

## Details

Points are located with a uniform bucket grid. Points outside the mesh
give `NA` (in `evaluate()`) or an error (in `evaluation_matrix()`).

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 4))
u <- interpolate(function(x, y) x^2, V)
evaluate(u, c(0.5, 0.25), c(0.5, 0.5), gradient = TRUE)
#>    value   dx dy
#> 1 0.2500 0.75  0
#> 2 0.0625 0.25  0
```
