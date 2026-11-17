# Solve a PDE problem

Generic solver entry point. Methods exist for elliptic, variational,
parabolic, reaction-diffusion and nonlinear problems.

## Usage

``` r
# S3 method for class 'elasticity_problem'
solve_pde(
  problem,
  method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
  ...
)

solve_pde(problem, ...)

# S3 method for class 'parabolic_problem'
solve_pde(
  problem,
  dt,
  method = c("crank-nicolson", "backward-euler"),
  output_times = NULL,
  solver = "auto",
  startup = 2L,
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

- method:

  For time-dependent problems: `"crank-nicolson"` (second order) or
  `"backward-euler"` (first order, L-stable).

- ...:

  Passed to methods: `method`, `tol` and other arguments of
  [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
  for linear problems, time stepping or Newton controls for the others.

- dt:

  Time step. The number of steps is `ceiling((T - t0) / dt)` and the
  step is then adjusted to land exactly on the final time.

- output_times:

  Times at which to store the solution; they must lie on the time grid.
  Default: every step.

- solver:

  Linear solver method, see
  [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md).

- startup:

  Number of initial backward Euler steps taken before Crank-Nicolson
  (Rannacher startup). They damp the oscillations that Crank-Nicolson
  produces for rough initial data or initial data that do not match the
  Dirichlet data; a fixed number of such steps keeps second order.
  Ignored for backward Euler.

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

## References

Rannacher, R. (1984). Finite element solution of diffusion problems with
irregular data. *Numerische Mathematik*, 43, 309-327.

## Examples

``` r
prob <- poisson_problem(pde_mesh_rectangle(nx = 8), source = 1, boundary = dirichlet_bc(0))
sol <- solve_pde(prob)
summary(sol)
#> Finite-element solution (Lagrange P1)
#>   mesh: 81 nodes, 128 triangles, h = 0.1768
#>   degrees of freedom: 81 (49 free)
#>   matrix nonzeros: 497
#>   solver: cholesky, status solved
#>   linear residual: 1.79e-16 (relative 1.64e-15)
#>   assembly 0.002 s, solve 0.001 s
#>   solution range [    0, 0.07278]
#>   (the linear residual measures the algebraic solve, not the discretization error)
```
