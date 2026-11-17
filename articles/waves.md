# Waves and elastodynamics

``` r

library(pdecompR)
```

## The damped wave equation

[`wave_problem()`](https://adamrsophie.github.io/pdecompR/reference/wave_problem.md)
poses \partial\_{tt} u + \gamma\\\partial_t u - \nabla\cdot(c^2\nabla u)
= f(x, y, t) with initial displacement and velocity and the usual
boundary conditions, whose data may depend on time. The
semidiscretization M\ddot{\mathbf u} + C\dot{\mathbf u} + K\mathbf u =
\mathbf F(t) is integrated by the Newmark average acceleration method
(\beta = 1/4, \gamma_N = 1/2). It is second order and unconditionally
stable. Without damping and forcing it conserves the discrete energy
\tfrac12\dot{\mathbf u}^\top M\dot{\mathbf u} + \tfrac12\mathbf u^\top
K\mathbf u exactly. The effective matrix is factorized once.

To isolate the time error, start from a discrete eigenmode K\mathbf x =
\lambda_h M\mathbf x (computed with
[`solve_eigen()`](https://adamrsophie.github.io/pdecompR/reference/solve_eigen.md)):
the semidiscrete solution is then exactly
\cos(\sqrt{\lambda_h}\\t)\\\mathbf x, and any remaining error comes from
time stepping.

``` r

mw <- pde_mesh_rectangle(nx = 16)
ev <- solve_eigen(elliptic_problem(mw, bcs = dirichlet_bc(0)), n = 1)
mode <- ev$vectors[[1]]; omega <- sqrt(ev$values[1])
p <- wave_problem(mw, initial = mode, bcs = dirichlet_bc(0), time = c(0, 1))
dts <- c(0.1, 0.05, 0.025, 0.0125)
err <- sapply(dts, function(dt) {
  s <- solve_pde(p, dt = dt, output_times = 1)
  max(abs(s$values[, 1] - cos(omega) * mode$values))
})
data.frame(dt = dts, error = err, rate = c(NA, log2(err[-4] / err[-1])))
#>       dt       error     rate
#> 1 0.1000 0.139091868       NA
#> 2 0.0500 0.035814070 1.957440
#> 3 0.0250 0.009019743 1.989368
#> 4 0.0125 0.002259093 1.997343
s <- solve_pde(p, dt = 0.05)
diff(range(s$diagnostics$energy)) / s$diagnostics$energy[1]
#> [1] 7.486982e-15
```

The relative energy drift is at rounding level.

Damping removes energy:

``` r

pd <- wave_problem(pde_mesh_rectangle(nx = 16), damping = 1,
                   initial = function(x, y) sin(pi * x) * sin(pi * y), bcs = dirichlet_bc(0), time = c(0, 3))
sd <- solve_pde(pd, dt = 0.05)
plot(seq(0, 3, by = 0.05), sd$diagnostics$energy, type = "l", xlab = "t", ylab = "energy")
```

![](waves_files/figure-html/unnamed-chunk-4-1.png)

## Elastodynamics

[`elastodynamics_problem()`](https://adamrsophie.github.io/pdecompR/reference/wave_problem.md)
adds inertia to an
[`elasticity_problem()`](https://adamrsophie.github.io/pdecompR/reference/elasticity_problem.md):
\rho\\\ddot{\mathbf u} + \gamma\rho\\\dot{\mathbf u} -
\nabla\cdot\sigma(\mathbf u) = \mathbf f. The static loads of the
elasticity problem are scaled by `load_factor(t)`. With damping, a
suddenly applied load leads to oscillations that settle at the static
solution:

``` r

m <- pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 16, ny = 4)
el <- elasticity_problem(m, E = 100, nu = 0.3, body_force = c(0, -1),
                         bcs = displacement_bc(c(0, 0), "left"))
dyn <- elastodynamics_problem(el, density = 1, damping = 1, time = c(0, 20))
sdyn <- solve_pde(dyn, dt = 0.1)
static <- solve_pde(el)
tip <- which.max(sdyn$space$dof_coords[, 1] + sdyn$space$dof_coords[, 2])
plot(sdyn$times, sdyn$values$uy[tip, ], type = "l", xlab = "t", ylab = "tip deflection")
abline(h = static$uy$values[tip], lty = 2)
```

![](waves_files/figure-html/unnamed-chunk-5-1.png)
