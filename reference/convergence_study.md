# Mesh convergence studies

Solves a problem on a sequence of meshes and reports errors and observed
convergence orders `log(e_{i-1} / e_i) / log(h_{i-1} / h_i)`.

## Usage

``` r
convergence_study(
  problem,
  exact,
  gradient = NULL,
  mesh_sizes = NULL,
  refinements = 3,
  norms = c("L2", "H1_semi"),
  ...
)
```

## Arguments

- problem:

  A stationary `pde_problem` (elliptic, nonlinear or elasticity).
  Variational problems are tied to their function space and are not
  supported.

- exact:

  Exact solution `function(x, y)`.

- gradient:

  Exact gradient (optional, see
  [`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)).

- mesh_sizes:

  Cell counts for rectangle meshes.

- refinements:

  Number of uniform refinements otherwise.

- norms:

  Norms to report.

- ...:

  Passed to
  [`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md).

## Value

A data frame of class `pde_convergence` with columns `h`, `ndof`, one
column per norm and the observed rates.
[`summary()`](https://rdrr.io/r/base/summary.html) reports the finest
errors and the last and mean observed rates per norm.

## Details

When the problem's mesh came from
[`pde_mesh_rectangle()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh_rectangle.md)
and `mesh_sizes` is given, meshes with `nx = mesh_sizes[i]` cells (and
the same aspect) are generated. Otherwise the problem's mesh is refined
uniformly `refinements` times.

Observed orders on coarse meshes can deviate from the asymptotic ones
(pre-asymptotic behavior). For P1 and a smooth solution the asymptotic
orders are 2 in L2 and 1 in H1; for P2 they are 3 and 2.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 4)
f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
prob <- poisson_problem(mesh, source = f, boundary = dirichlet_bc(0))
convergence_study(prob, function(x, y) sin(pi * x) * sin(pi * y),
                  mesh_sizes = c(4, 8, 16))
#> Mesh convergence study
#>        h ndof       L2 H1_semi rate_L2 rate_H1_semi
#>  0.35355   25 0.079077  0.8385      NA           NA
#>  0.17678   81 0.021133  0.4318   1.904       0.9575
#>  0.08839  289 0.005377  0.2175   1.974       0.9891
```
