# Nonlinear stationary problems

`nonlinear_problem()` defines the quasilinear problem

## Usage

``` r
nonlinear_problem(
  mesh,
  diffusion = 1,
  diffusion_derivative = NULL,
  reaction = NULL,
  reaction_derivative = NULL,
  source = 0,
  bcs = list(),
  degree = 1,
  initial = NULL,
  quad_degree = NULL
)

p_laplace_problem(
  mesh,
  p,
  source = 0,
  bcs = list(),
  epsilon = 1e-06,
  degree = 1,
  initial = NULL,
  quad_degree = NULL
)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- diffusion:

  `k(u)`: a number or `function(u)`.

- diffusion_derivative:

  `k'(u)` as `function(u)`, or `NULL`.

- reaction:

  `g(u)`: `NULL` or `function(u)`.

- reaction_derivative:

  `g'(u)` as `function(u)`, or `NULL`.

- source:

  Source `f`.

- bcs:

  A `pde_bc` or list of them. Boundary edges without a condition receive
  the natural condition `k du/dn = 0`.

- degree:

  Lagrange degree (1 or 2).

- initial:

  Optional initial guess (number, function or
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)).

- quad_degree:

  Quadrature exactness degree (default `2 degree + 2`).

- p:

  Exponent of the p-Laplacian (`p > 1`).

- epsilon:

  Regularization parameter of the p-Laplacian.

## Value

An object of class `nonlinear_problem`.

## Details

`-div(k(u) grad u) + g(u) = f`

and `p_laplace_problem()` the regularized p-Laplacian

`-div((|grad u|^2 + eps^2)^((p - 2) / 2) grad u) = f`.

Both are solved by Newton's method on the residual `F(u) = 0` of the
weak form, with the exact (consistent) Jacobian assembled from the
derivatives `k'(u)` and `g'(u)` (finite differences if not supplied),
and a backtracking line search on the residual norm. The convergence
history is returned; Newton is not assumed to converge, and failure is
reported.

The initial guess is the solution of the linear problem with `k(0)` and
`g'(0)` (for the p-Laplacian, the Laplace problem), unless `initial` is
given.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
prob <- nonlinear_problem(mesh, diffusion = function(u) 1 + u^2,
                          diffusion_derivative = function(u) 2 * u,
                          source = 10, bcs = dirichlet_bc(0))
sol <- solve_pde(prob)
sol$newton
#>   iteration     residual step
#> 1         0 2.011630e-01   NA
#> 2         1 8.790741e-03    1
#> 3         2 1.406079e-05    1
#> 4         3 2.957619e-11    1
#> 5         4 1.159439e-15    1
```
