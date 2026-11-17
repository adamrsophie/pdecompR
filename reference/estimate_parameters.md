# Parameter estimation for PDE models

Estimates parameters `theta` of a linear elliptic model from point
observations by minimizing

## Usage

``` r
estimate_parameters(
  problem_fun,
  theta0,
  observations,
  mesh = NULL,
  lower = -Inf,
  upper = Inf,
  regularization = 0,
  theta_ref = theta0,
  control = list(),
  derivative = NULL
)
```

## Arguments

- problem_fun:

  `function(theta)` or `function(theta, mesh)` returning an
  [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md).

- theta0:

  Initial parameter values.

- observations:

  A
  [`pde_observations()`](https://adamrsophie.github.io/pdecompR/reference/pde_observations.md)
  object.

- mesh:

  Optional mesh passed as second argument of `problem_fun`.

- lower, upper:

  Bounds on the parameters.

- regularization:

  Tikhonov weight `lambda`.

- theta_ref:

  Reference parameters of the regularization.

- control:

  Passed to [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- derivative:

  Optional exact operator derivatives, see
  [`pde_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/pde_sensitivity.md)
  and
  [`affine_derivative()`](https://adamrsophie.github.io/pdecompR/reference/affine_derivative.md);
  indices refer to the full parameter vector.

## Value

A list of class `pde_estimate` with the estimate, standard errors, the
`fixed` parameters, and `identifiability`: singular values, condition
number and right singular vectors of the weighted observation
sensitivity matrix of the free parameters, and their correlation matrix.

## Details

`1/2 sum_i ((u_h(x_i; theta) - y_i) / sd_i)^2 + lambda/2 ||theta - theta_ref||^2`

subject to the discretized PDE. The gradient is computed with the
discrete adjoint method
([`pde_sensitivity()`](https://adamrsophie.github.io/pdecompR/reference/pde_sensitivity.md));
the optimizer is [`stats::optim()`](https://rdrr.io/r/stats/optim.html)
with method `"L-BFGS-B"`, so bounds are supported. Forward solves are
cached so the objective and gradient at the same `theta` share one
solve.

At the estimate, the Gauss-Newton approximation `H = S^T W S + lambda I`
(with `S` the forward sensitivities of the predicted observations) gives
approximate standard errors. Without regularization they are
`sqrt(diag(H^-1))`; with `lambda > 0` they are the propagated
measurement noise `sqrt(diag(H^-1 S^T W S H^-1))` (the regularized
estimator is biased and this ignores the bias). They reflect measurement
noise only, under a local linearization; they do not include
discretization error. To expose the latter, the model is re-solved at
the estimate on a uniformly refined mesh (when `problem_fun` accepts a
mesh) and the change in the predicted observations is reported.

A local optimizer is used; for nonconvex problems the result is a local
minimum, and no global optimality is claimed.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
pf <- function(theta, mesh) poisson_problem(mesh, diffusion = theta[1], source = 1,
                                            boundary = dirichlet_bc(0))
truth <- solve_pde(pf(0.5, mesh))
xo <- c(0.3, 0.5, 0.7); yo <- c(0.5, 0.5, 0.4)
obs <- pde_observations(xo, yo, evaluate(truth$u, xo, yo), sd = 0.01)
est <- estimate_parameters(pf, theta0 = 1, observations = obs, mesh = mesh,
                           lower = 0.01)
est$estimate
#> [1] 0.5
```
