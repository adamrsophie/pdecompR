# Linear elasticity

Small-strain linear elasticity in two dimensions: find the displacement
`u = (u_x, u_y)` with

## Usage

``` r
elasticity_problem(
  mesh,
  E = 1,
  nu = 0.3,
  body_force = c(0, 0),
  bcs = list(),
  degree = 1,
  plane = c("strain", "stress"),
  quad_degree = NULL
)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- E, nu:

  Young's modulus and Poisson's ratio: numbers or `function(x, y)`.
  `E > 0` and `-1 < nu < 1/2` (plane strain) or `-1 < nu < 1` (plane
  stress).

- body_force:

  Body force `f`: a length-2 vector or a `function(x, y)` returning a
  two-column matrix.

- bcs:

  A list of
  [`displacement_bc()`](https://adamrsophie.github.io/pdecompR/reference/displacement_bc.md)
  and
  [`traction_bc()`](https://adamrsophie.github.io/pdecompR/reference/displacement_bc.md)
  conditions.

- degree:

  Polynomial degree of each displacement component (1 or 2).

- plane:

  `"strain"` or `"stress"`.

- quad_degree:

  Quadrature degree (default `2 * degree + 2`).

## Value

An object of class `elasticity_problem`;
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
returns a `pde_vector_solution` with the components `ux` and `uy`.

## Details

`-div sigma(u) = f`, `sigma(u) = lambda div(u) I + 2 mu eps(u)`,
`eps(u) = (grad u + grad u^T) / 2`,

with prescribed displacements on some boundary parts and prescribed
tractions `sigma(u) n = g` on others. Edges without a condition are
traction free. The displacement uses a vector-valued continuous Lagrange
space: each component lies in the scalar P1 or P2 space of the mesh, and
the global unknowns are ordered with all `x` components first.

The weak form is `int sigma(u) : eps(v) = int f . v + int_N g . v`; its
matrix is symmetric positive definite once the displacement conditions
remove the rigid motions (two translations and one rotation), so sparse
Cholesky applies. A problem without enough displacement conditions to
remove the rigid motions is rejected.

Lamé parameters follow from Young's modulus `E` and Poisson's ratio
`nu`: `mu = E / (2 (1 + nu))` and `lambda = E nu / ((1 + nu)(1 - 2 nu))`
for plane strain, or `lambda = E nu / (1 - nu^2)` for plane stress.
Low-order elements lock for nearly incompressible materials (`nu` close
to 1/2 in plane strain); P2 elements reduce but do not remove the
effect.

## Examples

``` r
mesh <- pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 16, ny = 4)
prob <- elasticity_problem(mesh, E = 100, nu = 0.3,
                           bcs = list(displacement_bc(c(0, 0), "left"),
                                      traction_bc(c(0, -1), "right")))
sol <- solve_pde(prob)
sol
#> <pde_vector_solution> P1 displacement, 170 dofs (160 free); solver cholesky (solved), relative residual 9.2e-13
#>   max |u| = 2.002
plot(sol)
```
