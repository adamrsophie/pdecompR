# Parameter sensitivities of linear elliptic problems

For a parameterized linear problem, the reduced discrete system after
Dirichlet elimination is `K(theta) u_f = r(theta)` with Dirichlet values
`u_d(theta)`. Differentiating gives the forward sensitivity equation

## Usage

``` r
pde_sensitivity(
  problem_fun,
  theta,
  mesh = NULL,
  objective = NULL,
  method = c("forward", "adjoint"),
  h = 1e-06,
  lower = -Inf,
  upper = Inf
)
```

## Arguments

- problem_fun:

  `function(theta)` or `function(theta, mesh)` returning an
  [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md).

- theta:

  Parameter vector.

- mesh:

  Optional mesh passed as second argument of `problem_fun`.

- objective:

  Optional `function(u)` of the full degree-of-freedom vector returning
  `list(value = , gradient = )`, for example from
  [`misfit_objective()`](https://adamrsophie.github.io/pdecompR/reference/pde_observations.md).

- method:

  `"forward"` (all solution sensitivities) or `"adjoint"` (objective
  gradient only).

- h:

  Relative difference step for the operator derivatives.

- lower, upper:

  Parameter bounds. Difference steps are shortened or made one-sided so
  that `problem_fun` is never called outside them.

## Value

A list of class `pde_sensitivity` with the `solution`, the
`sensitivities` matrix (forward), and the objective `value` and
`gradient` when `objective` is given.

## Details

`K du_f/dtheta_j = dr/dtheta_j - dK/dtheta_j u_f`,

valid when `K` is invertible and the data are differentiable. For a
scalar objective `J(u)` the discrete adjoint method solves
`K^T lambda = dJ/du_f` once and returns
`dJ/dtheta_j = dJ/du_d . du_d/dtheta_j + lambda . (dr/dtheta_j - dK/dtheta_j u_f)`,
at a cost independent of the number of parameters (apart from the
operator derivatives).

These are discrete sensitivities (differentiate, then the discrete
system is exact); they coincide with continuous adjoints only up to
discretization error.

The operator derivatives `dK/dtheta_j`, `dr/dtheta_j` and
`du_d/dtheta_j` are obtained by central differences of the assembled
operators with step `h_j = h max(1, |theta_j|)`. This is exact up to
rounding when the coefficients are affine in `theta`, and second-order
accurate otherwise. The solves themselves are not differenced.

## References

Hinze, M., Pinnau, R., Ulbrich, M. and Ulbrich, S. (2009). *Optimization
with PDE Constraints*. Springer.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
pf <- function(theta) poisson_problem(mesh, diffusion = theta[1], source = 1,
                                      boundary = dirichlet_bc(0))
s <- pde_sensitivity(pf, theta = 2)
# exact: u = u1 / theta, so du/dtheta = -u / theta
max(abs(s$sensitivities[, 1] + s$solution$u$values / 2))
#> [1] 1.496903e-12
```
