# Lagrange finite elements on triangles

`fe_element()` describes the continuous Lagrange element of degree 1
(P1) or 2 (P2) on the reference triangle `(0, 0)`, `(1, 0)`, `(0, 1)`.
With barycentric coordinates `l1 = 1 - xi - eta`, `l2 = xi`, `l3 = eta`:

## Usage

``` r
fe_element(degree = 1)
```

## Arguments

- degree:

  Polynomial degree, 1 or 2.

## Value

An object of class `fe_element` with fields `degree`, `ndof`, `nodes`
(reference coordinates of the local nodes) and functions
`basis(xi, eta)` and `grad(xi, eta)`. `basis` returns a matrix with one
row per point and one column per local basis function; `grad` returns a
list with the `xi` and `eta` derivatives in the same layout.

## Details

- P1: `phi_i = l_i`, `i = 1, 2, 3`, with nodes at the vertices.

- P2: `phi_i = l_i (2 l_i - 1)` at the vertices and `phi_4 = 4 l1 l2`,
  `phi_5 = 4 l2 l3`, `phi_6 = 4 l3 l1` at the midpoints of local edges
  1, 2 and 3.

Both satisfy the nodal property `phi_i(x_j) = delta_ij`.

## Examples

``` r
el <- fe_element(1)
el$basis(el$nodes[, 1], el$nodes[, 2])   # identity matrix
#>        xi eta
#> [1,] 1  0   0
#> [2,] 0  1   0
#> [3,] 0  0   1
```
