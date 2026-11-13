# Discretize a PDE problem

Builds the finite-element space and assembles the global system of a
linear stationary problem. The result separates the numerical
representation from the mathematical problem definition.

## Usage

``` r
# S3 method for class 'elasticity_problem'
discretize(problem, ...)

discretize(problem, ...)
```

## Arguments

- problem:

  A `pde_problem`.

- ...:

  Unused.

## Value

An object of class `pde_discretization` with the space, the full matrix
`A` and load `b` (before Dirichlet elimination), the Dirichlet degrees
of freedom and values, and assembly timings.

## Examples

``` r
prob <- poisson_problem(pde_mesh_rectangle(nx = 4), source = 1, boundary = dirichlet_bc(0))
disc <- discretize(prob)
disc
#> <pde_discretization> P1, 25 dofs (9 free), nnz(A) = 137
dim(disc$A)
#> [1] 25 25
```
