# Finite elements, spaces and quadrature

``` r

library(pdecompR)
```

## The reference triangle

All element computations happen on the reference triangle \hat T with
vertices (0,0), (1,0), (0,1). A physical triangle T with vertices x_1,
x_2, x_3 is its image under the affine map

F_T(\hat x) = x_1 + B_T \hat x, \qquad B_T = \[\\ x_2 - x_1 \\\\ x_3 -
x_1 \\\].

Counterclockwise orientation means \det B_T \> 0, and \|T\| = \det B_T /
2. Gradients transform with the inverse transpose, \nabla \phi =
B_T^{-T} \hat\nabla \hat\phi, and integrals pick up the factor \det B_T.

## Lagrange basis functions

With barycentric coordinates \lambda_1 = 1 - \xi - \eta, \lambda_2 =
\xi, \lambda_3 = \eta:

- **P1**: \phi_i = \lambda_i, nodes at the vertices;
- **P2**: \phi_i = \lambda_i(2\lambda_i - 1) at the vertices and
  4\lambda_1\lambda_2, 4\lambda_2\lambda_3, 4\lambda_3\lambda_1 at the
  midpoints of the local edges.

Both satisfy the nodal property \phi_i(x_j) = \delta\_{ij}:

``` r

p1 <- fe_element(1)
p1$basis(p1$nodes[, 1], p1$nodes[, 2])
#>        xi eta
#> [1,] 1  0   0
#> [2,] 0  1   0
#> [3,] 0  0   1
p2 <- fe_element(2)
round(p2$basis(p2$nodes[, 1], p2$nodes[, 2]), 12)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    0    0    0
#> [2,]    0    1    0    0    0    0
#> [3,]    0    0    1    0    0    0
#> [4,]    0    0    0    1    0    0
#> [5,]    0    0    0    0    1    0
#> [6,]    0    0    0    0    0    1
```

Gradients are available analytically and are checked against finite
differences in the test suite:

``` r

g <- p2$grad(0.2, 0.3)
rbind(xi = g$xi, eta = g$eta)
#>                z              
#> [1,] -1 -0.2 0.0  1.2 1.2 -1.2
#> [2,] -1  0.0 0.2 -0.8 0.8  0.8
```

The closed-form P1 element matrices serve as independent references:

K_T = \|T\|\\ G G^\top, \qquad M_T =
\frac{\|T\|}{12}\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix},

where the rows of G are the constant gradients of the three basis
functions.

``` r

v <- rbind(c(0, 0), c(2, 0), c(0, 1))
element_stiffness_p1(v)
#>       [,1]  [,2] [,3]
#> [1,]  1.25 -0.25   -1
#> [2,] -0.25  0.25    0
#> [3,] -1.00  0.00    1
element_mass_p1(v)
#>            [,1]       [,2]       [,3]
#> [1,] 0.16666667 0.08333333 0.08333333
#> [2,] 0.08333333 0.16666667 0.08333333
#> [3,] 0.08333333 0.08333333 0.16666667
```

## Function spaces and degrees of freedom

[`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md)
builds the global continuous space. P1 degrees of freedom are the mesh
nodes; P2 adds one per edge, numbered after the nodes.

``` r

mesh <- pde_mesh_rectangle(nx = 3)
V1 <- fe_space(mesh, 1); V2 <- fe_space(mesh, 2)
c(P1 = V1$ndof, P2 = V2$ndof, nodes = nrow(mesh$nodes), edges = nrow(mesh$edges))
#>    P1    P2 nodes edges 
#>    16    49    16    33
head(V2$cell_dofs, 3)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    2    6   17   35   26
#> [2,]    2    3    7   18   36   27
#> [3,]    3    4    8   19   37   28
boundary_dofs(V2, "bottom")
#> [1]  1  2  3  4 17 18 19
```

## Finite-element functions

[`interpolate()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
computes the nodal interpolant;
[`evaluate()`](https://adamrsophie.github.io/pdecompR/reference/evaluate.md)
returns values and gradients at arbitrary points (points are located
with a bucket grid).

``` r

u <- interpolate(function(x, y) exp(x) * cos(2 * y), V2)
evaluate(u, c(0.3, 0.7), c(0.2, 0.9), gradient = TRUE)
#>        value         dx        dy
#> 1  1.2433568  1.2598067 -1.056521
#> 2 -0.4522316 -0.4265419 -3.892735
c(exact = exp(0.3) * cos(0.4))
#>    exact 
#> 1.243302
plot(u)
```

![](finite-elements_files/figure-html/unnamed-chunk-7-1.png)

[`evaluation_matrix()`](https://adamrsophie.github.io/pdecompR/reference/evaluate.md)
returns the sparse matrix P with P\mathbf U equal to point values; it is
the observation operator of inverse problems.

## Quadrature

Element integrals use Gaussian rules on the reference triangle.
`quadrature_triangle(d)` is exact for polynomials of total degree d: 1
point for d \le 1, 3 points for d = 2, the 7-point Radon rule for d \le
5, and collapsed Gauss-Legendre products beyond.
`quadrature_interval(n)` is the n-point Gauss-Legendre rule, exact to
degree 2n - 1.

The exact monomial integral over \hat T is \int\_{\hat T} \xi^i \eta^j =
i!\\j!/(i+j+2)!:

``` r

q <- quadrature_triangle(5)
q
#> <fe_quadrature> triangle rule with 7 points, exact to degree 5
c(rule = sum(q$weights * q$points[, 1]^2 * q$points[, 2]^3),
  exact = factorial(2) * factorial(3) / factorial(7))
#>        rule       exact 
#> 0.002380952 0.002380952
```

The assembly default is degree 2p + 2 for elements of degree p: exact
for mass matrices with affine coefficients and for stiffness matrices
with quadratic coefficients. Error norms use degree 2p + 4. A rule that
is too coarse for the data would add a quadrature error to the
discretization error, so the defaults are deliberately generous.

[`integrate_mesh()`](https://adamrsophie.github.io/pdecompR/reference/integrate_mesh.md)
integrates any function over a mesh:

``` r

integrate_mesh(pde_mesh_rectangle(nx = 4), function(x, y) x^3 * y^2, degree = 5)
#> [1] 0.08333333
1 / 12
#> [1] 0.08333333
```
