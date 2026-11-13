# Finite-element functions

`fe_function()` wraps a coefficient vector in a space. `interpolate()`
returns the nodal (Lagrange) interpolant of a function.

## Usage

``` r
fe_function(space, values, name = "u")

interpolate(f, space, name = "u")
```

## Arguments

- space:

  An
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

- values:

  Numeric vector of length `space$ndof`.

- name:

  Optional label used in printing and plots.

- f:

  A number or vectorized function `function(x, y)`.

## Value

An object of class `fe_function`.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 4))
u <- interpolate(function(x, y) x + 2 * y, V)
evaluate(u, 0.3, 0.4)
#> [1] 1.1
```
