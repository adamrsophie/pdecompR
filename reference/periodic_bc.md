# Periodic boundary conditions

Identifies the degrees of freedom on boundary `to` with those on
boundary `from` translated by `offset`, so that `u(x + offset) = u(x)`
(and, in the weak form, the fluxes through the two boundaries cancel).
Pairs are matched by coordinates; every degree of freedom on `to` must
have a partner on `from`, so the two boundaries need matching nodes, as
in
[`pde_mesh_rectangle()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh_rectangle.md).
Several periodic conditions (for example in `x` and `y`) may be
combined, and corners are identified transitively. Supported for
[`elliptic_problem()`](https://adamrsophie.github.io/pdecompR/reference/elliptic_problem.md);
with no Dirichlet condition and no reaction the solution is fixed by a
zero-mean gauge.

## Usage

``` r
periodic_bc(from, to, offset)
```

## Arguments

- from, to:

  Boundary markers.

- offset:

  Translation `c(dx, dy)` taking `from` onto `to`.

## Value

An object of class `pde_bc`.

## Examples

``` r
m <- pde_mesh_rectangle(nx = 16)
p <- elliptic_problem(m, source = function(x, y) sin(2 * pi * x),
                      bcs = list(periodic_bc("left", "right", c(1, 0)),
                                 dirichlet_bc(0, c("bottom", "top"))))
range(solve_pde(p)$u$values)
#> [1] -0.02305718  0.02305718
```
