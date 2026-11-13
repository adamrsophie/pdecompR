# Adaptive finite-element solution

Runs the loop solve, estimate, mark, refine for a linear elliptic
problem, using
[`estimate_error()`](https://adamrsophie.github.io/pdecompR/reference/estimate_error.md),
[`mark_elements()`](https://adamrsophie.github.io/pdecompR/reference/mark_elements.md)
and newest vertex bisection
([`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)).

## Usage

``` r
adaptive_solve(
  problem,
  theta = 0.5,
  max_iterations = 8,
  max_dofs = 1e+05,
  tol = 0,
  exact = NULL,
  gradient = NULL,
  ...
)
```

## Arguments

- problem:

  An
  [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md).

- theta:

  Doerfler marking parameter.

- max_iterations:

  Maximum number of refinement steps.

- max_dofs:

  Stop when the number of degrees of freedom exceeds this.

- tol:

  Stop when the estimator falls below this.

- exact, gradient:

  Optional exact solution and gradient; when given, the true errors are
  recorded in the history.

- ...:

  Passed to
  [`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md).

## Value

A list of class `pde_adaptive` with the final `solution`, the final
`estimate` and a `history` data frame.
