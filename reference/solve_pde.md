# Solve a PDE problem

Generic solver entry point. Methods exist for elliptic, variational,
parabolic, reaction-diffusion and nonlinear problems.

## Usage

``` r
solve_pde(problem, ...)

# S3 method for class 'parabolic_problem'
solve_pde(
  problem,
  dt,
  method = c("crank-nicolson", "backward-euler"),
  output_times = NULL,
  solver = "auto",
  ...
)

# S3 method for class 'nonlinear_problem'
solve_pde(problem, tol = 1e-10, maxit = 30, method = "auto", ...)

# S3 method for class 'reaction_diffusion_problem'
solve_pde(
  problem,
  dt,
  method = c("imex-euler", "implicit-euler"),
  output_times = NULL,
  solver = "auto",
  newton_tol = 1e-10,
  newton_maxit = 25,
  ...
)
```

## Arguments

- problem:

  A `pde_problem`.

- ...:

  Passed to methods: `method`, `tol` and other arguments of
  [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
  for linear problems, time stepping or Newton controls for the others.

- dt:

  Time step. The number of steps is `ceiling((T - t0) / dt)` and the
  step is then adjusted to land exactly on the final time.

- method:

  For time-dependent problems: `"crank-nicolson"` (second order) or
  `"backward-euler"` (first order, L-stable).

- output_times:

  Times at which to store the solution; they must lie on the time grid.
  Default: every step.

- solver:

  Linear solver method, see
  [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md).

- tol, maxit:

  Newton tolerance on the free residual norm (relative to the initial
  residual, with an absolute floor of `1e-12`) and maximum number of
  Newton iterations.

- newton_tol, newton_maxit:

  Newton tolerance (on the residual norm relative to its initial value,
  or absolute `1e-12`) and maximum iterations per step.

## Value

A `pde_solution` (stationary problems) or `pde_transient`
(time-dependent problems).
