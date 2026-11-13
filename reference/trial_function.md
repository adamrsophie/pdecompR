# Weak-form language

A small language for linear scalar variational problems.
`trial_function()` and `test_function()` create the symbols `u` and `v`;
`grad()` and `dot()` build gradients and inner products; integrands are
combined with `*`, `+` and `-`, and turned into forms by
`integrate_domain()` (`dx`) and `integrate_boundary()` (`ds` over marked
boundary edges).

## Usage

``` r
trial_function(space)

test_function(space)

grad(x)

vector_coefficient(x)

dot(a, b)

integrate_domain(expr)

integrate_boundary(expr, marker = "boundary")
```

## Arguments

- space:

  An
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

- x, a, b:

  Weak-form expressions or coefficients.

- expr:

  A weak-form integrand.

- marker:

  Boundary markers for `integrate_boundary()`.

## Value

Expressions of class `fe_expr`, forms of class `fe_form`.

## Details

Supported integrands, with `c` a scalar coefficient (number or
`function(x, y)`) and `b` a vector coefficient (length-2 vector or
`function(x, y)` returning a two-column matrix):

- bilinear, domain: `c * dot(grad(u), grad(v))`, `c * u * v`,
  `dot(b, grad(u)) * v`, `u * dot(b, grad(v))`;

- bilinear, boundary: `c * u * v`;

- linear, domain: `c * v`, `dot(b, grad(v))`;

- linear, boundary: `c * v`.

Scalar coefficients multiply; products are evaluated pointwise at the
quadrature points. Write a coefficient function first and multiply it
afterwards (`f <- function(x, y) x * y; f * v`): in
`function(x, y) x * y * v` the `* v` belongs to the function body. The
functions `grad()` and `dot()` only act on these symbolic objects and do
not mask any base R function.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 8))
u <- trial_function(V); v <- test_function(V)
a <- integrate_domain(dot(grad(u), grad(v)))
f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
L <- integrate_domain(f * v)
prob <- variational_problem(a, L, dirichlet_bc(0))
sol <- solve_pde(prob)
```
