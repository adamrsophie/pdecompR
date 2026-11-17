# Adjoint gradients for nonlinear and time-dependent problems

`nonlinear_sensitivity()` returns the gradient of `J(u(theta))` for a
[`nonlinear_problem()`](https://adamrsophie.github.io/pdecompR/reference/nonlinear_problem.md)
family. At the converged state the discrete residual satisfies
`F(u, theta) = 0`; one solve with the transposed Newton Jacobian,
`J_F^T lambda = dJ/du`, gives
`dJ/dtheta = dJ/du_d du_d/dtheta - lambda' dF/dtheta`, where the
parameter derivative of the residual is a central difference at the
fixed state (no extra nonlinear solves).

## Usage

``` r
nonlinear_sensitivity(problem_fun, theta, objective, h = 1e-06, ...)

transient_sensitivity(
  problem_fun,
  theta,
  objective,
  dt,
  method = c("backward-euler", "crank-nicolson"),
  h = 1e-06,
  checkpoint = NULL,
  ...
)
```

## Arguments

- problem_fun:

  `function(theta)` returning a
  [`nonlinear_problem()`](https://adamrsophie.github.io/pdecompR/reference/nonlinear_problem.md)
  or
  [`parabolic_problem()`](https://adamrsophie.github.io/pdecompR/reference/parabolic_problem.md).

- theta:

  Parameter vector.

- objective:

  `function(u)` of the full degree-of-freedom vector (the final state
  for transient problems) returning `list(value = , gradient = )`, for
  example from
  [`misfit_objective()`](https://adamrsophie.github.io/pdecompR/reference/pde_observations.md).

- h:

  Relative difference step for parameter derivatives.

- ...:

  Passed to
  [`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md).

- dt, method:

  Time step and scheme for transient problems.

- checkpoint:

  For transient problems, `NULL` to store every state of the forward
  run, or a step interval `k`: only every `k`-th state is kept and the
  states in between are recomputed segment by segment during the
  backward sweep (one extra forward solve per step). With
  `k = ceiling(sqrt(N))` for `N` steps, memory falls from `N` to about
  `2 sqrt(N)` state vectors. The gradient is the same up to rounding.

## Value

A list with `value`, `gradient`, the `solution` and the adjoint.

## Details

`transient_sensitivity()` returns the gradient of an objective of the
final state of a
[`parabolic_problem()`](https://adamrsophie.github.io/pdecompR/reference/parabolic_problem.md)
family solved by backward Euler or Crank-Nicolson (without startup
steps). The forward run stores the states; a backward sweep with the
transposed step matrix computes the discrete adjoint, and the gradient
collects the derivatives of the operator, the mass matrix (which depends
on the parameters under SUPG stabilization), the loads and the initial
condition. It is the exact gradient of the discrete time-stepping scheme
(up to the difference steps of the operator derivatives), at the cost of
one extra solve per step, independent of the number of parameters.
Coefficients must not depend on time and the Dirichlet data must not
depend on `theta`. With `checkpoint`, only every `k`-th state is stored
and the others are recomputed during the backward sweep.

Check either gradient with Taylor remainders, as in
[`taylor_test()`](https://adamrsophie.github.io/pdecompR/reference/taylor_test.md).

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 6)
pf <- function(th) nonlinear_problem(mesh, diffusion = function(u) th[1] + u^2, source = th[2],
                                     bcs = dirichlet_bc(0))
ob <- misfit_objective(fe_space(mesh), pde_observations(0.5, 0.5, 0.05))
nonlinear_sensitivity(pf, c(1, 1), ob)$gradient
#> [1] -0.001576247  0.001578797
```
