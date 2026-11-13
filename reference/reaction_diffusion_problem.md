# Reaction-diffusion systems

`reaction_diffusion_problem()` defines a system of `S` species

## Usage

``` r
reaction_diffusion_problem(
  mesh,
  species = "u",
  diffusion = 1,
  reaction,
  jacobian = NULL,
  initial,
  bcs = list(),
  time = c(0, 1),
  degree = 1
)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- species:

  Character vector of species names.

- diffusion:

  Diffusion coefficients: a number, or a list or vector with one entry
  (number or `function(x, y)`) per species.

- reaction:

  `function(u, x, y, t)` where `u` is a matrix with one column per
  species (named) and one row per degree of freedom; it must return a
  matrix of the same shape.

- jacobian:

  Optional `function(u, x, y, t)` returning an array of dimension
  `c(n, S, S)` with entries `dR_a / du_b`.

- initial:

  A list with one initial condition per species.

- bcs:

  A named list of boundary conditions per species; species without an
  entry get zero-flux conditions.

- time:

  Start and end time.

- degree:

  Lagrange degree.

## Value

An object of class `reaction_diffusion_problem`.

## Details

`du_s/dt - div(D_s grad u_s) = R_s(u_1, ..., u_S, x, y, t)`, `s = 1..S`,

with per-species diffusion, initial conditions and boundary conditions.

The reaction term is discretized by the product (group finite-element)
approximation: `R` is evaluated at the degrees of freedom, interpolated,
and integrated against the test functions, giving `M R(U)`. This is
exact for linear reactions and second-order accurate in space for smooth
nonlinear ones.

Time integration (argument `method` of
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)):

- `"imex-euler"`: implicit diffusion, explicit reaction (first order);
  each species matrix is factored once.

- `"implicit-euler"`: fully implicit backward Euler with Newton's method
  on the coupled system and a backtracking line search. The reaction
  Jacobian is `jacobian` if given, otherwise computed by finite
  differences species by species (pointwise, so it is cheap and exact up
  to the difference error).

Positivity is not guaranteed by either scheme for general reactions.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 6)
prob <- reaction_diffusion_problem(
  mesh, species = "u", diffusion = 0.1,
  reaction = function(u, x, y, t) u * (1 - u),
  initial = list(u = function(x, y) 0.1 + 0.8 * x), time = c(0, 1))
sol <- solve_pde(prob, dt = 0.1)
```
