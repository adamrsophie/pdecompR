# Stresses and strains of an elasticity solution

Computes the strain `eps = (grad u + grad u^T) / 2`, the stress
`sigma = lambda tr(eps) I + 2 mu eps` and the von Mises stress
`sqrt(((s_xx - s_yy)^2 + (s_yy - s_zz)^2 + (s_zz - s_xx)^2) / 2 + 3 s_xy^2)`
of a solution of
[`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md).
The out-of-plane stress is `s_zz = nu (s_xx + s_yy)` in plane strain and
zero in plane stress.

## Usage

``` r
elasticity_stress(solution, at = c("elements", "nodes"))
```

## Arguments

- solution:

  A `pde_vector_solution` from
  [`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
  on an
  [`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md).

- at:

  `"elements"` or `"nodes"`.

## Value

For `"elements"`, a data frame with one row per element and columns `x`,
`y` (centroid), `exx`, `eyy`, `exy`, `sxx`, `syy`, `sxy`, `szz` and
`von_mises`. For `"nodes"`, a named list of P1
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
objects with the same fields.

## Details

Displacement gradients are discontinuous across elements. With
`at = "elements"` the fields are evaluated at element centroids (exact
for P1, where they are constant per element). With `at = "nodes"` the
centroid values are averaged over the elements around each vertex with
area weights and returned as P1 functions, which suits plotting and VTK
output; this recovery is not a higher-order stress approximation.

## Examples

``` r
mesh <- pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 16, ny = 4)
prob <- elasticity_problem(mesh, E = 100, nu = 0.3,
                           bcs = list(displacement_bc(0, "left"),
                                      traction_bc(c(0, -1), "right")))
st <- elasticity_stress(solve_pde(prob))
summary(st$von_mises)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.865   2.469   4.371   5.178   6.983  15.479 
```
