# Spatially varying diffusion tensors

Builds a symmetric diffusion tensor `K = [[xx, xy], [xy, yy]]` whose
entries are numbers or functions of `(x, y)` (or `(x, y, t)`), for
anisotropic media whose principal directions or strengths vary in space.
Use it wherever a diffusion coefficient is accepted. Positive
definiteness is checked at every quadrature point during assembly.

## Usage

``` r
tensor_coefficient(xx, xy = 0, yy = xx)
```

## Arguments

- xx, xy, yy:

  Tensor entries.

## Value

An object of class `pde_tensor`.

## Examples

``` r
# Diffusion aligned with circles around the origin
K <- tensor_coefficient(function(x, y) 1 + y^2, function(x, y) -x * y,
                        function(x, y) 1 + x^2)
p <- elliptic_problem(pde_mesh_rectangle(nx = 8), diffusion = K, source = 1,
                      bcs = dirichlet_bc(0))
max(solve_pde(p)$u$values)
#> [1] 0.05803095
```
