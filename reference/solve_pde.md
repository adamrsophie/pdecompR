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
  method = c("crank-nicolson", "backward-euler", "bdf2", "sdirk2"),
  output_times = NULL,
  solver = "auto",
  startup = 2L,
  adaptive = FALSE,
  tol = 1e-04,
  callback = NULL,
  restart = NULL,
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
  callback = NULL,
  restart = NULL,
  ...
)

# S3 method for class 'stokes_problem'
solve_pde(problem, method = c("lu", "gmres"), tol = 1e-10, ...)

# S3 method for class 'wave_problem'
solve_pde(
  problem,
  dt,
  output_times = NULL,
  solver = "auto",
  callback = NULL,
  restart = NULL,
  ...
)

# S3 method for class 'elastodynamics_problem'
solve_pde(
  problem,
  dt,
  output_times = NULL,
  solver = "auto",
  callback = NULL,
  restart = NULL,
  ...
)
```

## Arguments

- problem:

  A `pde_problem`.

- method:

  For time-dependent problems: `"crank-nicolson"` (second order),
  `"backward-euler"` (first order, L-stable), `"bdf2"` (second order,
  L-stable two-step method started with one backward Euler step) or
  `"sdirk2"` (Alexander's two-stage, second-order, L-stable singly
  diagonally implicit Runge-Kutta method). BDF2 and SDIRK2 require
  coefficients that do not depend on time (data may).

- ...:

  Passed to methods: `method`, `tol` and other arguments of
  [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
  for linear problems, time stepping or Newton controls for the others.

- dt:

  Time step. The number of steps is `ceiling((T - t0) / dt)` and the
  step is then adjusted to land exactly on the final time.

- output_times:

  Times at which to store the solution; they must lie on the time grid.
  Default: every step. With `adaptive = TRUE` they may be any times in
  the interval, and the steps land on them exactly.

- solver:

  Linear solver method, see
  [`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md).

- startup:

  Number of initial backward Euler steps taken before Crank-Nicolson
  (Rannacher startup). They damp the oscillations that Crank-Nicolson
  produces for rough initial data or initial data that do not match the
  Dirichlet data; a fixed number of such steps keeps second order.
  Ignored for backward Euler.

- adaptive:

  For `"sdirk2"`: choose the step size from an embedded first-order
  error estimate (L2 norm through the mass matrix), starting from `dt`.
  Steps whose estimate exceeds `tol` are rejected and retried.

- tol:

  Tolerance. Time-dependent problems: error tolerance per step for
  adaptive stepping, relative to `max(1, ||u||)`. Nonlinear problems:
  Newton tolerance on the free residual norm, relative to the initial
  residual with an absolute floor of `1e-12`. Stokes problems: GMRES
  relative residual tolerance.

- callback:

  Optional `function(t, u)` called after every accepted time step with
  the time and the full vector of nodal values (a matrix with one column
  per species for reaction-diffusion systems; the displacement vector
  for waves and elastodynamics). Use it to stream results to disk or to
  monitor long runs instead of storing every step.

- restart:

  A time-dependent solution from an earlier call. The integration
  continues from its final state (stored in `$state`: time, values and
  the extra state of multistep and Newmark schemes, the adaptive step
  size and the number of steps already taken) to `time[2]` of `problem`;
  the initial condition of `problem` is then ignored. With a fixed step
  size, a run split by restarts reproduces one uninterrupted run up to
  rounding. Adaptive SDIRK2 resumes with the stored next step size; it
  matches an uninterrupted run that also stops at the restart time
  (include it in `output_times`), and otherwise differs by about the
  error tolerance because the step sequences differ. Save the solution
  with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) for
  checkpointing. The state records its finite-element space (degree and
  degree-of-freedom coordinates in their numbering); a restart on a
  different mesh, degree or numbering is rejected. Reaction-diffusion
  states are matched to the problem's species by name. Loads and
  boundary data may change between runs.

- maxit:

  Maximum number of Newton iterations.

- newton_tol, newton_maxit:

  Newton tolerance (on the residual norm relative to its initial value,
  or `1e-14` relative to the size of the terms of the step equations)
  and maximum iterations per step.

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
#>   linear residual: 1.55e-16 (relative 1.42e-15)
#>   assembly 0.002 s, solve 0.002 s
#>   solution range [    0, 0.07278]
#>   (the linear residual measures the algebraic solve, not the discretization error)
```
