# Sparse linear solvers

`solve_linear_system()` solves `A x = b` for a sparse matrix and reports
diagnostics. `linear_solver()` prepares a reusable solver (for example a
Cholesky factorization) for repeated right-hand sides, as in time
stepping.

## Usage

``` r
solve_linear_system(
  A,
  b,
  method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
  tol = 1e-10,
  maxit = 10 * length(b),
  preconditioner = c("jacobi", "none"),
  restart = 50,
  x0 = NULL
)

linear_solver(
  A,
  method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
  tol = 1e-10,
  maxit = NULL,
  preconditioner = c("jacobi", "none"),
  restart = 50
)
```

## Arguments

- A:

  A square sparse matrix (`Matrix` class) or base matrix.

- b:

  Right-hand side vector.

- method:

  Solver, see Details.

- tol:

  Relative residual tolerance for iterative methods.

- maxit:

  Maximum iterations for iterative methods.

- preconditioner:

  `"jacobi"` or `"none"`.

- restart:

  Restart length for GMRES.

- x0:

  Initial guess for iterative methods.

## Value

`solve_linear_system()` returns a list with `x` and `info`, a list of
diagnostics (`method`, `status`, `residual`, `relative_residual`,
`iterations`, `tol`, `n`, `nnz`, `symmetric`, `time`). `linear_solver()`
returns a list with a function `solve(b)` and `info`.

## Details

Methods:

- `"auto"`: sparse Cholesky (CHOLMOD through `Matrix`) for symmetric
  matrices, falling back to sparse LU when the factorization shows the
  matrix is not positive definite; sparse LU for nonsymmetric matrices.

- `"cholesky"`: sparse Cholesky; fails for matrices that are not
  symmetric positive definite.

- `"lu"`: sparse LU (any nonsingular matrix).

- `"cg"`: preconditioned conjugate gradients (symmetric positive
  definite matrices only).

- `"bicgstab"`, `"gmres"`: preconditioned BiCGSTAB and restarted GMRES
  for general matrices.

Iterative methods use Jacobi (diagonal) preconditioning by default and
iterate on the scaled system `A y = b / max|b|`, which avoids overflow
and underflow. They report `"converged"` only when the true residual of
the returned solution satisfies `||b - A x|| <= tol * ||b||`, and
`"not converged"` otherwise. Direct methods report the same relative
residual, status `"inaccurate"` when the normwise backward error
`||b - A x|| / (||A|| ||x|| + ||b||)` exceeds `1e-8`, and status
`"near-singular"` when the ratio of the smallest to the largest pivot of
the factorization is below about `2e-13`. That ratio is a conditioning
diagnostic: the matrix is singular or nearly singular to working
precision, for example a stiffness matrix without boundary conditions. A
system with no unknowns returns method `"none"`. A small residual says
the linear system was solved accurately; it says nothing about the
discretization error of the PDE.

## References

Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, 2nd
edition. SIAM.

## Examples

``` r
V <- fe_space(pde_mesh_rectangle(nx = 8))
A <- assemble_matrix(V, "stiffness") + assemble_matrix(V, "mass")
b <- assemble_vector(V, 1)
s1 <- solve_linear_system(A, b)
s2 <- solve_linear_system(A, b, method = "cg", tol = 1e-12)
max(abs(s1$x - s2$x))
#> [1] 9.214851e-15
s2$info$iterations
#> [1] 23
```
