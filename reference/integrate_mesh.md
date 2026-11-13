# Integrate a function over a mesh

Applies a triangle quadrature rule element by element.

## Usage

``` r
integrate_mesh(mesh, f, degree = 4)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- f:

  A vectorized function `function(x, y)`.

- degree:

  Exactness degree of the quadrature rule.

## Value

The approximate integral.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 2)
integrate_mesh(mesh, function(x, y) x * y)   # exact: 1/4
#> [1] 0.25
```
