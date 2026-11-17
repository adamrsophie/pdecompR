# Taylor remainder test of an objective gradient

Checks a gradient `g` of `J(theta)` by the remainders
`|J(theta + e d) - J(theta) - e g . d|` for decreasing `e`. They must
decrease at second order (observed rates near 2) for a correct gradient
and only at first order for a wrong one.

## Usage

``` r
taylor_test(
  problem_fun,
  theta,
  objective,
  direction = NULL,
  mesh = NULL,
  eps = 10^-(1:4),
  method = c("adjoint", "forward"),
  derivative = NULL,
  lower = -Inf,
  upper = Inf
)
```

## Arguments

- problem_fun, mesh, objective, lower, upper, derivative:

  As in
  [`pde_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/pde_sensitivity.md).

- theta:

  Parameter vector.

- direction:

  Perturbation direction (default: all ones).

- eps:

  Step sizes, decreasing.

- method:

  Gradient method, `"adjoint"` or `"forward"`.

## Value

A data frame with `eps`, `remainder` and observed `rate`.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 6)
pf <- function(theta) poisson_problem(mesh, diffusion = theta[1], source = theta[2],
                                      boundary = dirichlet_bc(0))
V <- fe_space(mesh)
obs <- pde_observations(0.5, 0.5, 0.1)
taylor_test(pf, c(1, 1), objective = misfit_objective(V, obs))
#>     eps    remainder      rate
#> 1 1e-01 4.771183e-15        NA
#> 2 1e-02 4.774978e-16 0.9996547
#> 3 1e-03 4.809131e-17 0.9969048
#> 4 1e-04 3.958032e-18 1.0845873
```
