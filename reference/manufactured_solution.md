# Manufactured solutions with symbolic derivatives

`manufactured_solution()` takes an exact solution `u(x, y)` (or
`u(x, y, t)`) and coefficients written as strings or 'symcompR'
expressions, differentiates them symbolically with 'symcompR' and
returns vectorized R functions for the solution, its gradient and the
source term

## Usage

``` r
manufactured_solution(u, diffusion = 1, advection = NULL, reaction = 0)
```

## Arguments

- u:

  Exact solution, a string such as `"sin(pi*x)*sin(pi*y)"` or a
  'symcompR' expression in `x`, `y` and optionally `t`.

- diffusion:

  Scalar diffusion coefficient `k(x, y)` (string, number or expression).

- advection:

  Optional advection field: a length-2 list or character vector of
  components.

- reaction:

  Reaction coefficient `c(x, y)`.

## Value

A list of class `pde_manufactured` with functions `exact`, `gradient`
(returning a two-column matrix), `source`, `diffusion`, `advection` (or
`NULL`) and `reaction`, the symbolic expressions in `symbolic`, and
`time_dependent`. Time-dependent functions take `(x, y, t)`.

## Details

`f = du/dt - div(k grad u) + b . grad u + c u`,

where the time derivative is present only when `u` contains `t`. The
results plug directly into
[`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md),
[`heat_problem()`](https://adamrsophie.github.io/pdecompR/reference/parabolic_problem.md),
[`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)
and
[`convergence_study()`](https://adamrsophie.github.io/pdecompR/reference/convergence_study.md),
which removes hand-derived source terms as a source of error in
verification studies.

'symcompR' is an optional dependency; the function stops with a message
when it is not installed.

## Examples

``` r
if (requireNamespace("symcompR", quietly = TRUE)) {
  ms <- manufactured_solution("sin(pi*x)*sin(pi*y)", diffusion = "1 + x")
  prob <- elliptic_problem(pde_mesh_rectangle(nx = 8), diffusion = ms$diffusion,
                           source = ms$source, bcs = dirichlet_bc(0))
  errornorm(solve_pde(prob), ms$exact, "H1_semi", gradient = ms$gradient)
}
#> [1] 0.4318526
#> attr(,"norm")
#> [1] "H1_semi"
#> attr(,"gradient")
#> [1] "exact"
```
