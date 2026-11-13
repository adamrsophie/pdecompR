# Variational problems

`variational_problem()` pairs a bilinear form `a(u, v)` and a linear
form `L(v)` from the weak-form language with Dirichlet conditions
(Neumann and Robin terms are written directly as boundary integrals).
Solving it assembles both forms, eliminates the Dirichlet degrees of
freedom and solves the linear system.

## Usage

``` r
variational_problem(bilinear, linear, bcs = list())
```

## Arguments

- bilinear:

  An `fe_form` containing the trial and test functions.

- linear:

  An `fe_form` containing only the test function.

- bcs:

  Dirichlet conditions
  ([`dirichlet_bc()`](https://adamrsophie.github.io/pdecompR/reference/dirichlet_bc.md)).

## Value

An object of class `variational_problem`.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 8))
u <- trial_function(V); v <- test_function(V)
a <- integrate_domain(dot(grad(u), grad(v)))
L <- integrate_domain(1 * v)
sol <- solve_pde(variational_problem(a, L, dirichlet_bc(0)))
max(sol$u$values)
#> [1] 0.07278263
```
