# Point observations and data misfit

`pde_observations()` stores point measurements `y_i` at `(x_i, y_i)`
with standard deviations `sd_i`. `misfit_objective()` returns the
weighted least-squares misfit
`J(u) = 1/2 sum_i ((P u)_i - y_i)^2 / sd_i^2` and its gradient
`P^T W (P u - y)`, where `P` is the
[`evaluation_matrix()`](https://adamrsophie.github.io/pdecompR/reference/evaluate.md).

## Usage

``` r
pde_observations(x, y, values, sd = 1)

misfit_objective(space, observations)
```

## Arguments

- x, y:

  Observation coordinates.

- values:

  Observed values.

- sd:

  Measurement standard deviations (recycled).

- space:

  An
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

- observations:

  A `pde_observations` object.

## Value

`pde_observations()` returns a data frame of class `pde_observations`;
`misfit_objective()` returns a function of the degree-of-freedom vector.
