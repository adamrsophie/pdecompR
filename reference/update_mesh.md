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

## Examples

``` r
prob <- poisson_problem(pde_mesh_rectangle(nx = 2), source = 1, boundary = dirichlet_bc(0))
finer <- update_mesh(prob, pde_mesh_rectangle(nx = 8))
finer
#> <elliptic_problem> on a mesh with 128 triangles, P1
#>   dirichlet on boundary
```
