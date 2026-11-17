# Algebraic multigrid preconditioner

Builds a smoothed aggregation V-cycle from the matrix alone, so it
applies to P2 spaces, unstructured meshes and systems without a mesh
hierarchy, where
[`multigrid_preconditioner()`](https://adamrsophie.github.io/pdecompR/reference/multigrid_preconditioner.md)
does not. Unknowns are grouped into aggregates of strongly connected
neighbors (`|a_ij| >= strength * sqrt(a_ii a_jj)`); the tentative
prolongation is piecewise constant on the aggregates and, with
`smoothed = TRUE`, is improved by one damped Jacobi step. Coarse
operators are Galerkin products, smoothing is symmetric Gauss-Seidel and
the coarsest level (at most `coarse_size` unknowns) is solved directly.
It targets scalar symmetric positive definite systems such as diffusion
problems; use it with `method = "cg"`.

## Usage

``` r
amg_preconditioner(
  strength = 0.08,
  smoothed = TRUE,
  smoothing = 1L,
  coarse_size = 300L
)
```

## Arguments

- strength:

  Strength of connection threshold in `[0, 1)`.

- smoothed:

  Smooth the tentative prolongation.

- smoothing:

  Gauss-Seidel sweeps before and after the coarse correction.

- coarse_size:

  Largest system solved directly.

## Value

An object of class `pde_preconditioner`; its setup time is recorded in
the solver diagnostics as part of the solve time.

## References

Vanek, P., Mandel, J. and Brezina, M. (1996). Algebraic multigrid by
smoothed aggregation for second and fourth order elliptic problems.
*Computing*, 56, 179-196.

## Examples

``` r
prob <- poisson_problem(pde_mesh_rectangle(nx = 24), source = 1,
                        boundary = dirichlet_bc(0), degree = 2)
sol <- solve_pde(prob, method = "cg", preconditioner = amg_preconditioner())
sol$diagnostics$solver$iterations
#> [1] 22
```
