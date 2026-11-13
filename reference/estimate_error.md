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
`+ sum_{Neumann/Robin E in T} h_E ||g - alpha u_h - k grad u_h . n||_E^2`
`+ sum_{Dirichlet E in T} h_E ||d/ds (g_D - u_h)||_E^2`.

Edges with no condition are zero-flux Neumann edges. An edge under any
Dirichlet condition is a Dirichlet edge (as in the solver), and natural
conditions on the same edge add up. The Dirichlet term measures how well
the space represents the boundary data (tangential derivative by central
differences). The diffusion coefficient in the flux jumps is evaluated
on each side of an edge separately, so coefficients that jump across
element edges are handled; finite-element coefficients are allowed. For
P1 with constant diffusion the divergence term vanishes; variable scalar
diffusion contributes `grad k . grad u_h`, computed with central
differences of `k`. The estimator is reliable and efficient up to
constants for diffusion-dominated problems (Verfurth 2013); it is not
robust for strongly advection-dominated problems.

## References

Verfurth, R. (2013). *A Posteriori Error Estimation Techniques for
Finite Element Methods*. Oxford University Press.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
sol <- solve_pde(poisson_problem(mesh, source = 1, boundary = dirichlet_bc(0)))
est <- estimate_error(sol)
est
#> <pde_error_estimate> total 0.2269 over 128 elements (max 0.03453)
```
