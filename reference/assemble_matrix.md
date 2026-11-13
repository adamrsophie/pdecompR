# Finite-element assembly

`assemble_matrix()` assembles a global sparse matrix from element
contributions using sparse triplets, which are summed into a compressed
column `dgCMatrix`. Element integrals use triangle quadrature exact to
`quad_degree` (default `2 p + 2` for degree `p`), which is exact for
constant and affine coefficients.

## Usage

``` r
assemble_matrix(
  space,
  type = c("stiffness", "mass", "advection"),
  coef = NULL,
  t = NULL,
  quad_degree = NULL
)

assemble_vector(space, f, t = NULL, quad_degree = NULL)

assemble_boundary_matrix(
  space,
  alpha = 1,
  marker = "boundary",
  t = NULL,
  quad_degree = NULL
)

assemble_boundary_vector(
  space,
  g,
  marker = "boundary",
  t = NULL,
  quad_degree = NULL
)
```

## Arguments

- space:

  An
  [`fe_space()`](https://adamrsophie.github.io/pdecompR/reference/fe_space.md).

- type:

  Matrix type.

- coef:

  Coefficient, see Details.

- t:

  Optional time passed to time-dependent coefficients.

- quad_degree:

  Quadrature exactness degree.

- f, g, alpha:

  Scalar coefficients (numbers or `function(x, y)`, or
  `function(x, y, t)` together with `t`).

- marker:

  Boundary markers.

## Value

A sparse `dgCMatrix` or a numeric vector.

## Details

Matrix types (`i` is the test function, `j` the trial function):

- `"stiffness"`: `K_ij = int k grad(phi_j) . grad(phi_i)`; `coef` is the
  diffusion `k` (number, function or 2 by 2 matrix for a constant
  tensor).

- `"mass"`: `M_ij = int c phi_j phi_i`; `coef` defaults to 1.

- `"advection"`: `C_ij = int (b . grad(phi_j)) phi_i`; `coef` is the
  velocity `b` (length-2 vector or `function(x, y)` returning a two
  column matrix).

`assemble_vector()` assembles `F_i = int f phi_i`.
`assemble_boundary_matrix()` and `assemble_boundary_vector()` assemble
`int_Gamma alpha phi_j phi_i ds` and `int_Gamma g phi_i ds` over
boundary edges with the given markers, using Gauss-Legendre rules on
each edge.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 4))
K <- assemble_matrix(V, "stiffness")
M <- assemble_matrix(V, "mass")
sum(M)                                  # area of the domain
#> [1] 1
max(abs(K %*% rep(1, V$ndof)))          # constants are in the kernel
#> [1] 0
```
