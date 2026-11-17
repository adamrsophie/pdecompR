# Stokes flow with Taylor-Hood elements

`stokes_problem()` defines the incompressible Stokes equations

## Usage

``` r
stokes_problem(mesh, viscosity = 1, body_force = c(0, 0), bcs = list())

velocity_bc(value = c(0, 0), marker = "boundary", components = c("x", "y"))
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- viscosity:

  Kinematic viscosity `nu > 0` (number or `function(x, y)`).

- body_force:

  Body force `f` (length-2 vector or `function(x, y)` returning a
  two-column matrix).

- bcs:

  List of `velocity_bc()` and
  [`traction_bc()`](https://adamrsophie.github.io/pdecompR/reference/displacement_bc.md)
  conditions.

- value:

  Velocity `c(ux, uy)` or `function(x, y)` returning a two-column
  matrix.

- marker:

  Boundary markers.

- components:

  Velocity components to prescribe.

## Value

An object of class `stokes_problem`. Its solution has components
`velocity` (a vector solution with `ux` and `uy`, usable with
[`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)),
`pressure` (an
[`fe_function()`](https://adamrsophie.github.io/pdecompR/reference/fe_function.md))
and `diagnostics` (including the L2 norm of `div u_h`).

## Details

`-div(nu grad u) + grad p = f`, `div u = 0`

for velocity `u` and pressure `p`, discretized with the inf-sup stable
Taylor-Hood pair: continuous P2 velocity and P1 pressure. Velocity
conditions are given with `velocity_bc()` (any components, data may vary
in space);
[`traction_bc()`](https://adamrsophie.github.io/pdecompR/reference/displacement_bc.md)
prescribes the pseudo-traction `(nu grad u - p I) n`, and boundary parts
with no condition get the do-nothing outflow condition (zero
pseudo-traction). When the velocity is prescribed on the whole boundary,
the pressure is fixed by a zero-mean gauge and the prescribed velocity
must have zero net flux through the boundary.

[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
solves the saddle-point system by sparse LU (`method = "lu"`, the
default) or by right-preconditioned GMRES with the block-diagonal
preconditioner `diag(A, M_p / nu)` (`method = "gmres"`), where `A` is
the velocity Laplacian (factorized once) and `M_p` the pressure mass
matrix; its iteration count does not grow with mesh refinement.

## References

Taylor, C. and Hood, P. (1973). A numerical solution of the
Navier-Stokes equations using the finite element technique. *Computers
and Fluids*, 1, 73-100. Elman, H., Silvester, D. and Wathen, A. (2014).
*Finite Elements and Fast Iterative Solvers*, 2nd edition. Oxford
University Press.

## Examples

``` r
# Lid-driven cavity
m <- pde_mesh_rectangle(nx = 8)
p <- stokes_problem(m, bcs = list(velocity_bc(c(0, 0)), velocity_bc(c(1, 0), "top")))
s <- solve_pde(p)
range(s$velocity$ux$values)
#> [1] -0.1789024  1.0000000
```
