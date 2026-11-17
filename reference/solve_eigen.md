# Eigenvalue problems

Computes the eigenpairs `K x = lambda M x` of the discrete operator of
an
[`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md)
(diffusion and reaction, with Dirichlet conditions removing the
prescribed degrees of freedom; Neumann and Robin terms enter `K`) or an
[`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md)
(vibration modes, with `M` the mass matrix of unit density), nearest to
the shift `sigma`. Periodic conditions reduce both matrices as in
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md),
and the returned eigenvectors satisfy them.

## Usage

``` r
solve_eigen(problem, n = 6, sigma = NULL, tol = 1e-10, maxit = 500)
```

## Arguments

- problem:

  An
  [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md)
  without advection, or an
  [`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md).
  Sources, Neumann data and body forces are ignored; only the operator
  matters.

- n:

  Number of eigenpairs.

- sigma:

  Shift; the eigenvalues nearest to it are returned.

- tol:

  Relative residual tolerance.

- maxit:

  Maximum number of iterations.

## Value

A list of class `pde_eigen` with `values` (ascending distance to
`sigma`), `vectors` (a list of
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)s,
or of displacement pairs for elasticity, normalized to `x'Mx = 1`),
`residuals` and `iterations`.

## Details

The method is shift-and-invert subspace iteration with Rayleigh-Ritz
projection: one sparse factorization of `K - sigma M`, a block of
`n + 8` vectors solved with that factorization at each iteration, and a
dense generalized eigenproblem of the block size. It stops when every
requested pair has normwise relative residual
`||K x - lambda M x|| / ((||K|| + |lambda| ||M||) ||x||)` below `tol`.
Problems without Dirichlet conditions have a zero eigenvalue; the
default shift is then slightly negative so that `K - sigma M` is
nonsingular.

## Examples

``` r
p <- elliptic_problem(pde_mesh_rectangle(nx = 16), bcs = dirichlet_bc(0))
ev <- solve_eigen(p, n = 4)
ev$values / pi^2   # close to 2, 5, 5, 8
#> [1] 2.019310 5.082918 5.130183 8.305433
```
