# Wave equation and elastodynamics

`wave_problem()` defines the damped wave equation

## Usage

``` r
wave_problem(
  mesh,
  speed = 1,
  damping = 0,
  source = 0,
  initial = 0,
  velocity = 0,
  bcs = list(),
  time = c(0, 1),
  degree = 1
)

elastodynamics_problem(
  elasticity,
  density = 1,
  damping = 0,
  initial = c(0, 0),
  velocity = c(0, 0),
  time = c(0, 1),
  load_factor = function(t) 1
)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- speed:

  Wave speed `c` (number or `function(x, y)`).

- damping:

  Damping coefficient `gamma >= 0`.

- source:

  Source `f`, a number or `function(x, y)` or `function(x, y, t)`.

- initial, velocity:

  Initial displacement and velocity (numbers, functions of `(x, y)` or
  finite-element functions).

- bcs:

  Boundary conditions as for
  [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md).

- time:

  Start and end time.

- degree:

  Polynomial degree.

- elasticity:

  An
  [`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md)
  giving the material, loads and boundary conditions (displacement data
  are held fixed in time).

- density:

  Mass density `rho > 0` (a number).

- load_factor:

  `function(t)` scaling the body force and tractions.

## Value

An object of class `wave_problem` or `elastodynamics_problem`; solve it
with
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
and a time step `dt`.

## Details

`d2u/dt2 + gamma du/dt - div(c^2 grad u) = f(x, y, t)`

with initial displacement and velocity and Dirichlet, Neumann or Robin
conditions (data may depend on time). `elastodynamics_problem()` defines
`rho d2u/dt2 + gamma rho du/dt - div sigma(u) = f` for the displacement
of an
[`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md),
with the static loads of that problem multiplied by `load_factor(t)`.

Both semidiscretizations `M a + C v + K u = F(t)` are integrated by the
Newmark average acceleration method (`beta = 1/4`, `gamma = 1/2`):
second order, unconditionally stable and, without damping and forcing,
exactly conserving the discrete energy `1/2 v'Mv + 1/2 u'Ku`, which is
reported as `diagnostics$energy`. The effective matrix is factorized
once.

## References

Newmark, N. M. (1959). A method of computation for structural dynamics.
*Journal of the Engineering Mechanics Division*, 85, 67-94.

## Examples

``` r
u0 <- function(x, y) sin(pi * x) * sin(pi * y)
p <- wave_problem(pde_mesh_rectangle(nx = 8), initial = u0, bcs = dirichlet_bc(0),
                  time = c(0, 1))
s <- solve_pde(p, dt = 0.02)
range(s$diagnostics$energy)
#> [1] 2.435855 2.435855
```
