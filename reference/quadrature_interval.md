# Numerical quadrature rules

`quadrature_interval()` returns an `n`-point Gauss-Legendre rule on the
reference interval `[0, 1]`, computed with the Golub-Welsch algorithm.
It integrates polynomials of degree `2 n - 1` exactly.

## Usage

``` r
quadrature_interval(n)

quadrature_triangle(degree)
```

## Arguments

- n:

  Number of Gauss points.

- degree:

  Polynomial degree to integrate exactly.

## Value

A list with `points` (matrix), `weights` and `degree`, of class
`fe_quadrature`.

## Details

`quadrature_triangle()` returns a rule on the reference triangle with
vertices `(0, 0)`, `(1, 0)` and `(0, 1)` (area 1/2) that is exact for
all polynomials of total degree `degree`:

- degree 0 or 1: the centroid rule (1 point);

- degree 2: the three-point interior rule (Strang and Fix);

- degree 3 to 5: the seven-point rule of Radon (Stroud T2:5-1);

- higher degrees: a collapsed (Duffy) tensor product of Gauss-Legendre
  rules with positive weights and interior points.

Rules are mapped to a physical triangle by the affine map
`x = x1 + B (xi, eta)`; weights are multiplied by `|det B|`.

## References

Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
rules. *Mathematics of Computation*, 23, 221-230.

Stroud, A. H. (1971). *Approximate Calculation of Multiple Integrals*.
Prentice-Hall.

## Examples

``` r
q <- quadrature_triangle(4)
sum(q$weights)                               # area of the reference triangle
#> [1] 0.5
sum(q$weights * q$points[, 1]^2 * q$points[, 2]^2)  # exact: 1/180
#> [1] 0.005555556
```
