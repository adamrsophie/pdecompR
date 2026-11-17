# Piecewise coefficients by region

`region_coefficient()` builds a coefficient that takes one value per
mesh region
([`mesh_regions()`](https://adamrsophie.github.io/pdecompR/reference/mesh_regions.md)),
for example the conductivities of different materials. The result is a
vectorized `function(x, y)`, so it can be used wherever a coefficient is
accepted (diffusion, reaction, source, Young's modulus, error norms). On
the mesh it was built from, assembly reads the element values directly;
on other meshes (refined or not) the value at a point is that of the
region containing it in the original mesh, so the coefficient keeps its
geometry under refinement.

## Usage

``` r
region_coefficient(mesh, values)
```

## Arguments

- mesh:

  The
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md)
  whose regions define the pieces.

- values:

  A named numeric vector or list with one value per region (all regions
  of the mesh must be covered).

## Value

A function of class `pde_region_coefficient`.

## Details

At points on an interface between regions the value of one adjacent
region is returned; quadrature points and the error estimator never
evaluate exactly on element edges.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
mesh <- set_region(mesh, "right", function(x, y) x > 0.5)
k <- region_coefficient(mesh, c(domain = 1, right = 100))
prob <- poisson_problem(mesh, diffusion = k, source = 1, boundary = dirichlet_bc(0))
range(solve_pde(prob)$u$values)
#> [1] 0.00000000 0.02859314
```
