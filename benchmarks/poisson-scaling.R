# Assembly and solve times for the P1 and P2 Poisson problem on the unit
# square. Run from the package root with
#   Rscript benchmarks/poisson-scaling.R
# Results depend on the machine and the BLAS; record both with the numbers.
library(pdecompR)
f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
run <- function(n, degree) {
  mesh <- pde_mesh_rectangle(nx = n)
  t_mesh <- system.time(mesh <- pde_mesh_rectangle(nx = n))[["elapsed"]]
  s <- solve_pde(poisson_problem(mesh, source = f, boundary = dirichlet_bc(0), degree = degree))
  data.frame(degree = degree, n = n, ndof = s$space$ndof, elements = nrow(mesh$triangles),
             mesh_s = t_mesh, assembly_s = s$diagnostics$assembly_time,
             solve_s = s$diagnostics$solve_time,
             nnz = s$diagnostics$solver$nnz,
             rel_residual = s$diagnostics$solver$relative_residual)
}
res <- rbind(do.call(rbind, lapply(c(32, 64, 128, 256), run, degree = 1)),
             do.call(rbind, lapply(c(32, 64, 128), run, degree = 2)))
print(res, digits = 3)
cat("\nR", R.version.string, "| Matrix", as.character(packageVersion("Matrix")),
    "|", Sys.info()[["sysname"]], Sys.info()[["machine"]], "\n")
