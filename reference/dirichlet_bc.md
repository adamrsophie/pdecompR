# Boundary conditions

Boundary conditions are selected by boundary marker, never by coordinate
comparisons.

## Usage

``` r
dirichlet_bc(value = 0, marker = "boundary")

neumann_bc(flux = 0, marker = "boundary")

robin_bc(alpha, value = 0, marker = "boundary")
```

## Arguments

- value:

  Boundary value.

- marker:

  Boundary markers; `"boundary"` selects the whole boundary.

- flux:

  Outward flux `k du/dn`.

- alpha:

  Robin coefficient (non-negative).

## Value

An object of class `pde_bc`.

## Details

- `dirichlet_bc(value, marker)`: `u = value` on the marked boundary.
  Enforced strongly by eliminating the boundary degrees of freedom and
  moving their contribution to the right-hand side, which keeps
  symmetric positive definite systems symmetric positive definite.

- `neumann_bc(flux, marker)`: `k du/dn = flux` (outward normal). Adds
  `int flux v ds` to the load.

- `robin_bc(alpha, value, marker)`: `k du/dn + alpha u = value`. Adds
  `int alpha u v ds` to the bilinear form and `int value v ds` to the
  load.

`value`, `flux` and `alpha` may be numbers, `function(x, y)` or, for
time-dependent problems, `function(x, y, t)`.

Dirichlet data are imposed by nodal interpolation at the boundary
degrees of freedom. When several Dirichlet conditions share a degree of
freedom (for example a corner), the condition listed last wins. A degree
of freedom with a Dirichlet condition ignores Neumann and Robin
contributions.

## Examples

``` r
dirichlet_bc(0, "boundary")
#> <pde_bc> dirichlet on boundary 
neumann_bc(function(x, y) x, "top")
#> <pde_bc> neumann on top 
robin_bc(alpha = 2, value = 1, marker = "right")
#> <pde_bc> robin on right 
```
