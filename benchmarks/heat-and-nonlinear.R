# Time-stepping and Newton benchmarks. Run from the package root with
#   Rscript benchmarks/heat-and-nonlinear.R
library(pdecompR)
u0 <- function(x, y) sin(pi * x) * sin(pi * y)
for (n in c(32, 64)) {
  p <- heat_problem(pde_mesh_rectangle(nx = n), initial = u0, bcs = dirichlet_bc(0), time = c(0, 0.1))
  s <- solve_pde(p, dt = 0.001, output_times = 0.1)
  cat(sprintf("heat P1 n = %3d dofs = %6d steps = %d  assembly %.2fs  stepping %.2fs\n", n, s$space$ndof,
              s$diagnostics$steps, s$diagnostics$assembly_time, s$diagnostics$solve_time))
}
k <- function(u) 1 + u^2
for (n in c(32, 64)) {
  p <- nonlinear_problem(pde_mesh_rectangle(nx = n), diffusion = k, diffusion_derivative = function(u) 2 * u,
                         source = 10, bcs = dirichlet_bc(0))
  t <- system.time(s <- solve_pde(p))[["elapsed"]]
  cat(sprintf("quasilinear n = %3d dofs = %6d Newton iterations %d  total %.2fs\n", n, s$space$ndof,
              nrow(s$newton) - 1L, t))
}
