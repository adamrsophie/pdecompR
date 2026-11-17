# Write a time series for ParaView

Writes every stored time of a time-dependent solution as an XML VTK file
(`.vtu`) and a ParaView collection file (`.pvd`) that lists them with
their times, so ParaView opens the series as one animated data set.

## Usage

``` r
write_pvd(x, file, name = "u")
```

## Arguments

- x:

  A `pde_transient` from
  [`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md).

- file:

  Path of the `.pvd` file; the `.vtu` files are written next to it with
  the same stem and the time index appended.

- name:

  Name of the point data array for single-species solutions.

## Value

The `.pvd` path, invisibly.

## Examples

``` r
p <- heat_problem(pde_mesh_rectangle(nx = 4), initial = function(x, y) x * y,
                  bcs = dirichlet_bc(0), time = c(0, 0.1))
s <- solve_pde(p, dt = 0.05)
f <- tempfile(fileext = ".pvd")
write_pvd(s, f)
```
