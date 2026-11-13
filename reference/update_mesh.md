# Replace the mesh of a problem

Returns a copy of `problem` posed on `mesh`; coefficients and boundary
conditions are kept. Used by convergence studies and adaptive loops.

## Usage

``` r
update_mesh(problem, mesh)
```

## Arguments

- problem:

  A `pde_problem`.

- mesh:

  The new
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

## Value

The updated problem.
