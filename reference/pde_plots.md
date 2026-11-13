# Plot meshes and finite-element functions

Base-graphics plots that work in any R session (the solvers themselves
never need a graphics device).

## Usage

``` r
# S3 method for class 'pde_mesh'
plot(x, markers = TRUE, main = "Mesh", ...)

# S3 method for class 'fe_function'
plot(
  x,
  type = c("fill", "contour"),
  col = grDevices::hcl.colors(64, "viridis"),
  nlevels = 12,
  mesh = FALSE,
  main = x$name,
  ...
)

# S3 method for class 'pde_solution'
plot(x, ...)

# S3 method for class 'pde_transient'
plot(x, time = max(x$times), species = 1, main = NULL, ...)

# S3 method for class 'pde_convergence'
plot(x, main = "Mesh convergence", ...)

plot_indicators(
  mesh,
  values,
  col = grDevices::hcl.colors(64, "inferno", rev = TRUE),
  main = "Element indicators",
  ...
)

# S3 method for class 'pde_adaptive'
plot(x, main = "Adaptive refinement", ...)
```

## Arguments

- x:

  Object to plot.

- markers:

  Color boundary edges by marker.

- main:

  Title.

- ...:

  Passed to the underlying plotting functions.

- type:

  `"fill"` or `"contour"`.

- col:

  Color palette (vector of colors).

- nlevels:

  Number of contour levels or color levels.

- mesh:

  Draw mesh edges on top of a filled plot.

- time:

  Time to plot for transient solutions.

- species:

  Species for reaction-diffusion solutions.

- values:

  Element-wise values for `plot_indicators()`.

## Value

The input, invisibly.

## Details

- `plot.pde_mesh()` draws the triangles and colors boundary edges by
  marker.

- `plot.fe_function()` draws a filled plot (each triangle colored by the
  mean of its vertex values) or contour lines from a regular-grid
  evaluation (`type = "contour"`).

- `plot.pde_solution()` plots the solution; `plot.pde_transient()` plots
  a stored time; `plot.pde_convergence()` draws errors against `h` on
  log scales with reference slopes; `plot_indicators()` draws
  element-wise values such as error indicators.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 6)
plot(mesh)

u <- interpolate(function(x, y) sin(pi * x) * y, fe_space(mesh))
plot(u)

plot(u, type = "contour")
```
