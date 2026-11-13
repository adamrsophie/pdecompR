# Mesh size

Mesh size

## Usage

``` r
mesh_size(mesh)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

## Value

The largest edge length `h` of the mesh.

## Examples

``` r
mesh_size(pde_mesh_rectangle(nx = 4))
#> [1] 0.3535534
```
