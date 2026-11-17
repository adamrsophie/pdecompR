# Errors against a reference solution

`errornorm()` computes the error between a finite-element function and
an exact solution by element-wise quadrature of a rule of higher degree
than the one used for assembly (default `2 p + 4` for degree `p`). Nodal
errors are not used: they measure a different (superconvergent)
quantity.

## Usage

``` r
errornorm(
  u,
  exact,
  type = c("L2", "H1_semi", "H1", "energy"),
  gradient = NULL,
  diffusion = 1,
  quad_degree = NULL,
  elementwise = FALSE
)
```

## Arguments

- u:

  An
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md),
  `pde_solution` or `pde_vector_solution` (from
  [`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md);
  then `exact` returns two columns and `gradient` is a list of two
  gradient functions, one per component).

- exact:

  Exact solution `function(x, y)`, or an
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
  reference solution on any mesh of the same domain, possibly of another
  degree; its gradient is then used for gradient norms. The integral
  runs on a partition on which both functions are polynomials: the finer
  of two nested meshes (for example from
  [`refine_uniform()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
  or
  [`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)),
  or the overlay of all pairwise triangle intersections for meshes that
  are not nested (slower). The result is exact up to quadrature and
  symmetric in the two functions.

- type:

  Norm type.

- gradient:

  Exact gradient `function(x, y)` returning a two-column matrix.

- diffusion:

  Diffusion coefficient for the energy norm: a number, a function of
  `(x, y)` or a constant 2 by 2 matrix.

- quad_degree:

  Quadrature exactness degree.

- elementwise:

  Return element contributions (squared) as attribute `elements`.

## Value

The error (a number) with attributes.

## Details

Norms:

- `"L2"`: `||u - u_h||_{L2}`;

- `"H1_semi"`: `|u - u_h|_{H1} = ||grad(u - u_h)||_{L2}`;

- `"H1"`: the full norm `sqrt(L2^2 + H1_semi^2)`;

- `"energy"`: `sqrt(int k |grad(u - u_h)|^2)` with diffusion `k`.

Gradient-based norms need the exact gradient. When `gradient` is not
given it is approximated by central differences with step `1e-5` times
the domain diameter; the result is then labeled with attribute
`gradient = "finite difference"`. Its accuracy depends on the smoothness
of `exact` (the truncation error is of order step squared times the
third derivative), and it is wrong near jumps of the gradient, so pass
`gradient` for verification work. The step does not change when the
domain is translated; on a domain far from the origin it is raised where
the coordinates cannot represent it, and the rounding error then grows
like `1e-16 |x|` divided by the step.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
sol <- solve_pde(poisson_problem(mesh, source = f, boundary = dirichlet_bc(0)))
errornorm(sol, function(x, y) sin(pi * x) * sin(pi * y), "L2")
#> [1] 0.02113281
#> attr(,"norm")
#> [1] "L2"
#> attr(,"gradient")
#> [1] "exact"
```
