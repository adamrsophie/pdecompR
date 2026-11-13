# Linear elliptic problems

`elliptic_problem()` defines the scalar second-order problem

## Usage

``` r
elliptic_problem(
  mesh,
  diffusion = 1,
  reaction = 0,
  advection = NULL,
  source = 0,
  bcs = list(),
  degree = 1,
  stabilization = c("none", "supg"),
  quad_degree = NULL
)

poisson_problem(
  mesh,
  diffusion = 1,
  source = 0,
  boundary = list(),
  degree = 1,
  bcs = boundary
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

- bcs:

  A `pde_bc` or list of them. Boundary edges without a condition receive
  the natural condition `k du/dn = 0`.

- degree:

  Lagrange degree (1 or 2).

- stabilization:

  `"none"` or `"supg"`.

- quad_degree:

  Quadrature exactness degree (default `2 degree + 2`).

- boundary:

  Alias of `bcs` for `poisson_problem()`.

## Value

An object of class `pde_problem` (subclass `elliptic_problem`).

## Details

`-div(k grad u) + b . grad u + c u = f` in the domain,

with Dirichlet, Neumann and Robin boundary conditions.
`poisson_problem()` is the special case `b = 0`, `c = 0`.

The weak form is: find `u` with the Dirichlet data such that for all
test functions `v` vanishing on the Dirichlet boundary

`int k grad u . grad v + (b . grad u) v + c u v + sum_Robin int alpha u v`
`= int f v + sum_Neumann int g v + sum_Robin int r v`.

Coefficients may be numbers or vectorized `function(x, y)`; `diffusion`
may also be a constant symmetric positive definite 2 by 2 matrix, and
`advection` a length-2 vector or a function returning a two-column
matrix. Uniform ellipticity (`k > 0`) is checked at the quadrature
points.

With `stabilization = "supg"` (P1 only) the streamline upwind
Petrov-Galerkin terms
`sum_T tau_T int_T (b . grad u + c u - f) (b . grad v)` are added, with
`tau_T = h_T / (2 |b|) (coth(Pe) - 1 / Pe)`, `Pe = |b| h_T / (2 k)`,
evaluated at the element centroid with `h_T` the longest edge. The
method is consistent: the exact solution satisfies the stabilized
equations.

Problems without Dirichlet or Robin conditions and without reaction are
pure Neumann problems: the solution is unique only up to a constant. The
solver then checks the compatibility condition `int f + int g = 0` and
returns the solution with zero mean, labeled as such. The discrete check
is `|sum_i b_i| <= compatibility_tol * sum_i |b_i|` (an argument of
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md),
default `1e-6`), which tolerates quadrature error in compatible data; a
Lagrange multiplier for the zero-mean constraint absorbs the remaining
defect and is reported.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
prob <- poisson_problem(mesh, source = function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y),
                        boundary = dirichlet_bc(0))
sol <- solve_pde(prob)
sol
#> <pde_solution> P1, 81 dofs; solver cholesky (solved), relative residual 2.06e-15
#>   u range [    0, 0.9872]
```
