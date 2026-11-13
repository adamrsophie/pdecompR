# Boundary conditions for linear elasticity

`displacement_bc()` prescribes displacement components on marked edges;
`components = "x"` fixes only `u_x` (a roller), for example.
`traction_bc()` prescribes the surface force `sigma(u) n`.

## Usage

``` r
displacement_bc(value = c(0, 0), marker = "boundary", components = c("x", "y"))

traction_bc(value = c(0, 0), marker = "boundary")
```

## Arguments

- value:

  A length-2 vector or a `function(x, y)` returning a two-column matrix.

- marker:

  Boundary markers; `"boundary"` selects the whole boundary.

- components:

  Displacement components to prescribe (`"x"`, `"y"`).

## Value

An object of class `pde_vector_bc`.

## Examples

``` r
displacement_bc(c(0, 0), "left")
#> <pde_vector_bc> displacement [x, y] = (    0,     0) on left
displacement_bc(0, "bottom", components = "y")
#> <pde_vector_bc> displacement [y] = (    0,     0) on bottom
traction_bc(function(x, y) cbind(0, -x), "top")
#> <pde_vector_bc> traction = function(x, y) on top
```
