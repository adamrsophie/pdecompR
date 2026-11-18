# Sparse assembly and linear solvers

``` r

library(pdecompR)
library(Matrix)
```

## Element loops without loops

Assembly computes, for every element T and every pair of local basis
functions,

A^T\_{ij} = \sum_q w_q \det B_T \left\[ c\\\phi_i\phi_j +
\phi_i\\(b\cdot\nabla\phi_j) + \nabla\phi_i^\top D\\\nabla\phi_j
\right\](\hat x_q),

with i the test index (row) and j the trial index (column). pdecompR
vectorizes over elements: the loops run over quadrature points and the
3\times3 (or 6\times6) local indices, while each operation acts on a
vector of length m (the number of elements). The local matrices are then
summed into a global sparse matrix through triplets (I, J, A^T\_{ij}),
which
[`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
accumulates into the compressed column format.

``` r

V <- fe_space(pde_mesh_rectangle(nx = 16))
K <- assemble_matrix(V, "stiffness")
M <- assemble_matrix(V, "mass")
C <- assemble_matrix(V, "advection", c(1, 0))
class(K); dim(K); nnzero(K)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
#> [1] 289 289
#> [1] 1377
```

The structure reflects the mesh connectivity:

``` r

image(K[1:60, 1:60], main = "stiffness pattern (first 60 dofs)")
```

![](assembly_files/figure-html/unnamed-chunk-4-1.png)

## Properties worth checking

``` r

max(abs(K %*% rep(1, V$ndof)))   # constants are in the kernel of K
#> [1] 0
sum(M)                           # 1' M 1 is the area
#> [1] 1
isSymmetric(K); isSymmetric(M)
#> [1] TRUE
#> [1] TRUE
range(eigen(as.matrix(M), only.values = TRUE)$values)
#> [1] 0.0002663858 0.0038627258
```

The stiffness matrix with Dirichlet rows and columns removed is
symmetric positive definite; the advection matrix is not symmetric.

## Tensor coefficients

Anisotropic diffusion uses a symmetric tensor
`tensor_coefficient(xx, xy, yy)`, whose entries are numbers or functions
of `(x, y)`. Positive definiteness is checked at every quadrature point.

``` r

D <- tensor_coefficient(xx = 1, xy = 0.3, yy = function(x, y) 0.2 + x)
D
#> <pde_tensor> xx =     1, xy =   0.3, yy = function
sa <- solve_pde(elliptic_problem(pde_mesh_rectangle(nx = 16), diffusion = D, source = 1,
                                 bcs = dirichlet_bc(0)))
range(sa$u$values)
#> [1] 0.00000000 0.08892322
```

## Load vectors and boundary terms

``` r

b <- assemble_vector(V, function(x, y) x * y)
sum(b)                                            # int x y = 1/4
#> [1] 0.25
sum(assemble_boundary_vector(V, 1, "left"))       # length of the left edge
#> [1] 1
```

## Projection between meshes

`l2_project(f, space)` computes the L^2 projection: the function in the
space with the same integrals against every basis function. Unlike
interpolation it preserves the integral of `f`. For a finite-element
function on an unrelated mesh, the load integrals are computed exactly
on the overlay (common refinement) of the two meshes, so the integral
over the whole domain is preserved up to rounding. Integrals over single
elements or groups of elements are not preserved in general; that would
need a discontinuous space. The two meshes must cover the same domain,
otherwise the projection is rejected.
[`errornorm()`](https://adamrsophie.github.io/pdecompR/reference/errornorm.md)
uses the same overlay to compare solutions on unrelated meshes.

``` r

a <- interpolate(function(x, y) sin(5 * x) * y, fe_space(pde_mesh_rectangle(nx = 12)))
W <- fe_space(pde_mesh_rectangle(nx = 5, pattern = "left"))
pr <- l2_project(a, W); it <- interpolate(a, W)
c(source = sum(assemble_vector(a$space, 1) * a$values),
  projection = sum(assemble_vector(W, 1) * pr$values),
  interpolation = sum(assemble_vector(W, 1) * it$values))
#>        source    projection interpolation 
#>    0.07003947    0.07003947    0.06626694
```

## Invariance under renumbering

Assembly does not depend on the order of nodes: renumbering the nodes
permutes the matrix.

``` r

m <- pde_mesh_rectangle(nx = 4)
perm <- sample(nrow(m$nodes)); inv <- order(perm)
m2 <- pde_mesh(m$nodes[perm, ], matrix(inv[m$triangles], ncol = 3))
K1 <- as.matrix(assemble_matrix(fe_space(m), "stiffness"))
K2 <- as.matrix(assemble_matrix(fe_space(m2), "stiffness"))
max(abs(K2 - K1[perm, perm]))
#> [1] 0
```

## Linear solvers

[`solve_linear_system()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
chooses sparse Cholesky for symmetric matrices (with a fallback to
sparse LU when the matrix is not positive definite) and sparse LU
otherwise. Preconditioned conjugate gradients, BiCGSTAB and restarted
GMRES are also available.

``` r

A <- K + M
rhs <- assemble_vector(V, 1)
for (m in c("auto", "lu", "cg", "bicgstab", "gmres")) {
  s <- solve_linear_system(A, rhs, method = m, tol = 1e-12)
  cat(sprintf("%-8s %-9s iterations %4s  relative residual %.2e\n", m, s$info$method,
              s$info$iterations, s$info$relative_residual))
}
#> auto     cholesky  iterations   NA  relative residual 1.82e-13
#> lu       lu        iterations   NA  relative residual 1.77e-13
#> cg       cg        iterations   54  relative residual 4.51e-13
#> bicgstab bicgstab  iterations   43  relative residual 5.04e-13
#> gmres    gmres     iterations   54  relative residual 8.73e-13
```

Every solve reports a status, the residual and its relative size, the
iteration count, the matrix size and the number of nonzeros. Iterative
methods report `"not converged"` instead of returning silently:

``` r

solve_linear_system(A, rhs, method = "cg", maxit = 3)$info$status
#> [1] "not converged"
```

For repeated solves with the same matrix (time stepping),
[`linear_solver()`](https://adamrsophie.github.io/pdecompR/reference/solve_linear_system.md)
factors once and returns a reusable
[`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html) function.
A matrix right-hand side solves several systems with one factorization.

Direct solves first scale the matrix symmetrically by D^{-1/2} A
D^{-1/2} with D = \|\mathrm{diag}(A)\|, which removes the effect of
badly scaled coefficients on the pivoting and on the residual test.

### Preconditioners and multigrid

The Krylov methods take `preconditioner = "jacobi"` (the default),
`"ssor"`, `"none"`, or a function `function(r)` that applies an
approximate inverse. On a hierarchy of nested meshes,
[`multigrid_preconditioner()`](https://adamrsophie.github.io/pdecompR/reference/multigrid_preconditioner.md)
builds a geometric multigrid V-cycle: exact P1 prolongation from the
refinement record, Galerkin coarse operators, symmetric Gauss-Seidel
smoothing and a direct solve on the coarsest mesh. Its iteration count
stays bounded under refinement, while Jacobi’s grows like h^{-1}:

``` r

meshes <- list(pde_mesh_rectangle(nx = 4))
for (k in 1:4) meshes[[k + 1]] <- refine_uniform(meshes[[k]])
p <- poisson_problem(meshes[[5]], source = 1, boundary = dirichlet_bc(0))
its <- sapply(list(jacobi = "jacobi", ssor = "ssor", multigrid = multigrid_preconditioner(meshes)),
              function(pc) solve_pde(p, method = "cg", preconditioner = pc)$diagnostics$solver$iterations)
its
#>    jacobi      ssor multigrid 
#>       131       104         9
```

The multigrid preconditioner works for P1 problems on meshes produced by
[`refine_uniform()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md)
or
[`refine_mesh()`](https://adamrsophie.github.io/pdecompR/reference/refine_uniform.md).

Without a mesh hierarchy, or for P2 spaces,
[`amg_preconditioner()`](https://adamrsophie.github.io/pdecompR/reference/amg_preconditioner.md)
builds the hierarchy from the matrix by smoothed aggregation. On an
unstructured L-shaped mesh with P2 elements its iteration count barely
moves while Jacobi’s doubles with each halving of `h` (at `h = 0.025`,
about 130,000 unknowns, Jacobi CG reaches its iteration limit while AMG
needs 37 iterations):

``` r

L <- pde_domain("polygon", vertices = cbind(c(0, 2, 2, 1, 1, 0), c(0, 0, 1, 1, 2, 2)))
sapply(c(0.1, 0.05), function(h) {
  p <- poisson_problem(pde_mesh_polygon(L, h = h), source = 1, boundary = dirichlet_bc(0), degree = 2)
  c(h = h, jacobi = solve_pde(p, method = "cg")$diagnostics$solver$iterations,
    amg = solve_pde(p, method = "cg", preconditioner = amg_preconditioner())$diagnostics$solver$iterations)
})
#>         [,1]   [,2]
#> h        0.1   0.05
#> jacobi 337.0 682.00
#> amg     32.0  35.00
```

In two dimensions the sparse Cholesky factorization remains faster than
AMG preconditioned CG at every size measured (0.53 s against 0.82 s for
66,000 P2 unknowns on the unit square, with similar peak memory), so the
direct solver stays the default. AMG is the choice among the iterative
methods: far fewer iterations than Jacobi or SSOR, and it applies where
geometric multigrid cannot (P2 spaces, unstructured meshes without a
refinement hierarchy).

## Cost

Assembly is linear in the number of elements, and so is the memory of
the sparse matrices. The supernodal Cholesky factorization of
two-dimensional finite-element matrices costs roughly O(N^{1.5})
operations with fill reducing orderings, which keeps meshes with 10^5 to
10^6 unknowns practical on a workstation.

``` r

for (n in c(16, 32, 64)) {
  mesh <- pde_mesh_rectangle(nx = n)
  t <- system.time(s <- solve_pde(poisson_problem(mesh, source = 1, boundary = dirichlet_bc(0))))
  cat(sprintf("n = %3d  dofs = %6d  assembly %.2fs  solve %.2fs\n", n, s$space$ndof,
              s$diagnostics$assembly_time, s$diagnostics$solve_time))
}
#> n =  16  dofs =    289  assembly 0.00s  solve 0.00s
#> n =  32  dofs =   1089  assembly 0.00s  solve 0.00s
#> n =  64  dofs =   4225  assembly 0.01s  solve 0.01s
```

Where coefficients change every time step or Newton iteration, the
matrix is reassembled on a fixed sparsity pattern: the map from element
entries to the nonzeros is computed once, and each reassembly is one
sparse matrix-vector product, about twelve times faster than building
the matrix from triplets. The element integrals themselves then dominate
(about 80 percent of a reassembly for P1 on a 128 by 128 mesh), so the
end-to-end gain is modest: 7 percent for a heat problem with
time-dependent diffusion and 4 percent for a quasilinear Newton solve.
Building the plan costs about as much as 9 direct assemblies for P1 and
20 for P2, so the first fills are assembled directly and the plan is
built only once that many fills have been requested: short runs do not
pay for it. Load vectors are reassembled in time only when the source or
the natural boundary data depend on time; moving Dirichlet values alone
do not trigger it.

`benchmarks/protocol.R` in the source repository times the clean public
workload (mesh creation and one
[`solve_pde()`](https://adamrsophie.github.io/pdecompR/reference/solve_pde.md)
call) in a fresh R process under `/usr/bin/time`, so the peak memory
belongs to that workload alone, times the assembly components separately
in another process, and computes the error of each case in a third
process, so the error quadrature enters neither the solve time nor its
peak memory. Each case runs twice to warm up and then five times; every
repeat is saved, the median and range are reported, and the commit,
BLAS, CPU and thread settings are recorded. P1, P2 and AMG solves are
compared at a fixed error. On the unit square P2 reaches an L2 error of
`1e-5` with 4,225 unknowns, where P1 needs more than 66,000.
