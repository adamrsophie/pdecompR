# A posteriori error estimation

`estimate_error()` computes the residual-based error estimator for a
linear elliptic solution:

## Usage

``` r
estimate_error(solution)
```

## Arguments

- solution:

  A `pde_solution` of an
  [`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md).

## Value

A list of class `pde_error_estimate` with `eta` (per element), `total`,
and the interior and edge contributions.

## Details

`eta_T^2 = h_T^2 ||f - c u_h - b . grad u_h + div(k grad u_h)||_T^2`
`+ 1/2 sum_{interior E in T} h_E ||[k grad u_h . n]||_E^2`
`+ sum_{Neumann/Robin E in T} h_E ||g - alpha u_h - k grad u_h . n||_E^2`.

Edges with no condition are zero-flux Neumann edges; Dirichlet edges do
not contribute. For P1 with constant diffusion the divergence term
vanishes; variable scalar diffusion contributes `grad k . grad u_h`,
computed with central differences of `k`. The estimator is reliable and
efficient up to constants for diffusion-dominated problems (Verfurth
2013); it is not robust for strongly advection-dominated problems.

## References

Verfurth, R. (2013). *A Posteriori Error Estimation Techniques for
Finite Element Methods*. Oxford University Press.
