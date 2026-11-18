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
evaluated at the element centroid with `h_T` the longest edge (for a
matrix `K`, `k` is the streamline diffusivity `b'Kb / |b|^2`). For P1
elements the diffusive part of the element residual vanishes, so the
method is consistent when the diffusion is constant; with variable
diffusion the term `grad k . grad u_h` is omitted, which keeps the
convergence order but is not strictly consistent.

A connected part of the domain (mesh components joined by periodic
conditions) on which no Dirichlet value, reaction term or Robin term
with positive `alpha` acts keeps the constants in its kernel: there the
solution is unique only up to a constant. The solver gauges each such
component separately, checks its compatibility condition (for a
symmetric operator `int f + int g = 0` over the component) and returns
the solution with zero mean on it, labeled as such. Components that are
anchored are solved as usual. The compatibility defect is read from the
Lagrange multiplier of the zero-mean constraint, after the load is
normalized, and compared with `compatibility_tol * sum_i |b_i|` over the
component (an argument of
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md),
default `1e-6`, at most `0.01`). For a nonsymmetric operator (advection)
the left null vector is not constant, and data that are compatible for
the PDE leave a discrete defect of the size of the discretization error.
A relative defect up to `(h / diameter)^2` is then accepted silently,
one up to one percent with a warning, and anything larger is an error on
every mesh: the discrete equations have no solution, and refining the
mesh is the remedy when the data are compatible. Whenever the defect is
above rounding level (`1e-10`), the returned function solves the
equations with the source shifted by the multiplier times the constant
mass, not the original ones: `diagnostics$gauge$modified` is then
`TRUE`, `diagnostics$gauge$original_problem_solved` is `FALSE`, and
`diagnostics$solver$relative_residual` is the residual of the original
equations.

## References

Brooks, A. N. and Hughes, T. J. R. (1982). Streamline
upwind/Petrov-Galerkin formulations for convection dominated flows.
*Computer Methods in Applied Mechanics and Engineering*, 32, 199-259.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
prob <- poisson_problem(mesh, source = function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y),
                        boundary = dirichlet_bc(0))
sol <- solve_pde(prob)
sol
#> <pde_solution> P1, 81 dofs; solver cholesky (solved), relative residual 1.54e-15
#>   u range [    0, 0.9872]
```
