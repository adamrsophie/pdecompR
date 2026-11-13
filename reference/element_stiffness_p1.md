# Local P1 element matrices

Closed-form P1 stiffness and mass matrices of a single triangle, used as
independent references for the quadrature-based assembly.

## Usage

``` r
element_stiffness_p1(vertices, k = 1)

element_mass_p1(vertices)
```

## Arguments

- vertices:

  A 3 by 2 matrix of vertex coordinates (counterclockwise).

- k:

  Constant diffusion coefficient.

## Value

A 3 by 3 matrix.

## Details

The stiffness matrix is `K_e = |T| G G^T` where the rows of `G` are the
constant gradients of the three basis functions; the mass matrix is
`M_e = |T| / 12 * (1 + I)` (2 on the diagonal, 1 elsewhere).

## Examples

``` r
v <- rbind(c(0, 0), c(1, 0), c(0, 1))
element_stiffness_p1(v)
#>      [,1] [,2] [,3]
#> [1,]  1.0 -0.5 -0.5
#> [2,] -0.5  0.5  0.0
#> [3,] -0.5  0.0  0.5
element_mass_p1(v)
#>            [,1]       [,2]       [,3]
#> [1,] 0.08333333 0.04166667 0.04166667
#> [2,] 0.04166667 0.08333333 0.04166667
#> [3,] 0.04166667 0.04166667 0.08333333
```
