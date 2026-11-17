# Geometric multigrid preconditioner

Builds a V-cycle preconditioner for P1 systems on a hierarchy of nested
meshes, each obtained from the previous one by one call of
[`refine_uniform()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
or
[`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md).
Prolongation is exact P1 interpolation (new nodes are edge midpoints,
taken from the refinement record), coarse operators are Galerkin
products `P^T A P`, smoothing is forward Gauss-Seidel before and
backward Gauss-Seidel after the coarse correction (so the cycle is
symmetric for symmetric `A` and can precondition CG), and the coarsest
level is solved directly. With Dirichlet conditions the fine-level
prolongation is restricted to the free degrees of freedom automatically.

## Usage

``` r
multigrid_preconditioner(meshes, smoothing = 2L)
```

## Arguments

- meshes:

  List of nested meshes, coarse to fine; the last one is the mesh of the
  system.

- smoothing:

  Gauss-Seidel sweeps before and after the coarse correction.

## Value

An object of class `pde_preconditioner`.

## Details

Pass it as `preconditioner` to
[`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md),
[`linear_solver()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
or
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
for elliptic problems, together with an iterative `method`. The
iteration count then stays nearly constant as the mesh is refined.

## Examples

``` r
m0 <- pde_mesh_rectangle(nx = 4)
ms <- list(m0)
for (k in 1:3) ms[[k + 1]] <- refine_uniform(ms[[k]])
prob <- poisson_problem(ms[[4]], source = 1, boundary = dirichlet_bc(0))
sol <- solve_pde(prob, method = "cg", preconditioner = multigrid_preconditioner(ms))
sol$diagnostics$solver$iterations
#> [1] 8
```
