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
