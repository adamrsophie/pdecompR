# Extract a time snapshot

Extract a time snapshot

## Usage

``` r
solution_at(solution, time, species = 1)
```

## Arguments

- solution:

  A `pde_transient` solution.

- time:

  A stored output time (the nearest stored time is used).

- species:

  For reaction-diffusion systems, the species name or index.

## Value

An
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md).

## Examples

``` r
prob <- heat_problem(pde_mesh_rectangle(nx = 4), initial = function(x, y) sin(pi * x) * sin(pi * y),
                     bcs = dirichlet_bc(0), time = c(0, 0.1))
sol <- solve_pde(prob, dt = 0.05)
u <- solution_at(sol, 0.1)
evaluate(u, 0.5, 0.5)
#> [1] 0.2175204
```
