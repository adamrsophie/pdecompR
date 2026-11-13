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
