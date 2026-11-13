# Distributed optimal control

Solves the linear-quadratic optimal control problem

## Usage

``` r
optimal_control(
  mesh,
  target,
  alpha,
  diffusion = 1,
  source = 0,
  bcs = dirichlet_bc(0),
  lower = NULL,
  upper = NULL,
  control = list()
)
```

## Arguments

- mesh:

  A
  [`pde_mesh()`](https://adamrsophie.github.io/pdecompR/reference/pde_mesh.md).

- target:

  Desired state `u_target(x, y)`.

- alpha:

  Control cost (positive).

- diffusion:

  Diffusion coefficient.

- source:

  Fixed source `f`.

- bcs:

  Boundary conditions of the state equation.

- lower, upper:

  Optional control bounds (numbers).

- control:

  Passed to [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list of class `pde_control` with the `control` and `state` functions,
the `objective`, the gradient norm and solver details.

## Details

`min_q 1/2 ||u - u_target||^2 + alpha/2 ||q||^2`

subject to `-div(k grad u) = f + q` with the given boundary conditions,
where the control `q` is a P1 function and norms are `L2` norms
(computed with the mass matrix). The problem is strictly convex for
`alpha > 0`.

Without bounds the first-order optimality (KKT) system is solved
directly with sparse LU. With bounds `lower <= q <= upper`, the reduced
problem is solved with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) (`"L-BFGS-B"`)
using adjoint gradients and a single Cholesky factorization of the state
operator.

## Examples

``` r
mesh <- pde_mesh_rectangle(nx = 8)
oc <- optimal_control(mesh, target = function(x, y) sin(pi * x) * sin(pi * y),
                      alpha = 1e-4)
oc
#> Distributed optimal control (KKT (sparse LU))
#>   objective 0.00479426, projected gradient norm 3.76e-17
#>   control range [-2.025e-11, 19.68]
```
