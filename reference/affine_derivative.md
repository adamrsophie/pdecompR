# Exact derivatives of an affine problem family

For problems whose assembled system is affine in the parameters,
`A(theta) = A0 + sum(theta_j A_j)` and
`b(theta) = b0 + sum(theta_j b_j)` (for example a diffusion, reaction or
source that is a linear combination of fixed fields), the components are
the unit differences `A_j = A(theta + e_j) - A(theta)`.
`affine_derivative()` assembles them once and returns a `derivative`
function for
[`pde_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/pde_sensitivity.md)
and
[`taylor_test()`](https://adamrsophie.github.io/pdecompR/reference/taylor_test.md);
later calls cost no assembly. The result is exact only for affine
families: verify with
[`taylor_test()`](https://adamrsophie.github.io/pdecompR/reference/taylor_test.md).

## Usage

``` r
affine_derivative(problem_fun, theta, mesh = NULL)
```

## Arguments

- problem_fun, mesh:

  As in
  [`pde_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/pde_sensitivity.md).

- theta:

  Parameter vector at which the components are assembled.

## Value

A `function(theta, j)` returning `list(A = , b = , dirichlet = )`.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 6)
pf <- function(theta) poisson_problem(mesh, diffusion = theta[1], source = theta[2],
                                      boundary = dirichlet_bc(0))
dv <- affine_derivative(pf, c(1, 1))
obs <- pde_observations(0.5, 0.5, 0.1)
pde_sensitivity(pf, c(2, 3), objective = misfit_objective(fe_space(mesh), obs),
                method = "adjoint", derivative = dv)$gradient
#> [1] -0.0004420534  0.0002947023
```
