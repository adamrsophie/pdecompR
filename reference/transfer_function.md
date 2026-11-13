# Transfer a finite-element function to another mesh

Interpolates `u` into the Lagrange space of the given degree on `mesh`
by evaluating `u` at the new degrees of freedom. When `mesh` was
produced from `u`'s mesh by
[`refine_uniform()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
or
[`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
and `u` is P1, the transfer is exact (new nodes are edge midpoints).

## Usage

``` r
transfer_function(u, mesh, degree = u$space$degree)
```

## Arguments

- u:

  An
  [`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md).

- mesh:

  Target
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- degree:

  Degree of the target space.

## Value

An
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md)
on `mesh`.
