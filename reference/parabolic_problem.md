# Linear parabolic problems

`parabolic_problem()` defines

## Usage

``` r
parabolic_problem(
  mesh,
  diffusion = 1,
  reaction = 0,
  advection = NULL,
  source = 0,
  initial = 0,
  bcs = list(),
  time = c(0, 1),
  degree = 1,
  stabilization = c("none", "supg"),
  quad_degree = NULL
)

heat_problem(
  mesh,
  diffusivity = 1,
  source = 0,
  initial = 0,
  bcs = list(),
  time = c(0, 1),
  degree = 1
)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- diffusion:

  Diffusion coefficient `k`.

- reaction:

  Reaction coefficient `c` (non-negative for coercivity).

- advection:

  Velocity `b`, or `NULL`.

- source:

  Source `f`.

- initial:

  Initial condition: number, `function(x, y)` or
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md).

- bcs:

  A `pde_bc` or list of them. Boundary edges without a condition receive
  the natural condition `k du/dn = 0`.

- time:

  Start and end time.

- degree:

  Lagrange degree (1 or 2).

- stabilization:

  `"none"` or `"supg"`.

- quad_degree:

  Quadrature exactness degree (default `2 degree + 2`).

- diffusivity:

  Diffusivity `kappa` of the heat equation.

## Value

An object of class `parabolic_problem`.

## Details

`du/dt - div(k grad u) + b . grad u + c u = f(x, y, t)` for `t` in
`(time[1], time[2])`, with `u(x, y, time[1]) = initial(x, y)` and
boundary conditions whose data may depend on time (`function(x, y, t)`).
The diffusion, advection and reaction coefficients are time independent,
so the system matrices are assembled once.

`heat_problem()` is the heat equation `du/dt - kappa Laplace(u) = f`.

The finite-element semidiscretization is `M du/dt + A u = F(t)` with
mass matrix `M` and operator matrix `A` (including Robin terms). Time
stepping uses the theta method; see
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
arguments `dt` and `method`.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
prob <- heat_problem(mesh, diffusivity = 1,
                     initial = function(x, y) sin(pi * x) * sin(pi * y),
                     bcs = dirichlet_bc(0), time = c(0, 0.1))
sol <- solve_pde(prob, dt = 0.01)
sol
#> <pde_transient> P1, 81 dofs, 10 steps of crank-nicolson (dt =  0.01)
#>   stored times: 11 in [    0,   0.1]
```
