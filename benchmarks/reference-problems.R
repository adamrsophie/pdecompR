# Reference benchmark suite (PROMPT.md section 27).
#
# Runs the ten minimum benchmark problems on a sequence of meshes and records
# mesh size, degrees of freedom, assembly and solve time, peak memory, error
# against the reference solution, observed convergence order, solver
# configuration and tolerances. Results go to
# benchmarks/results/reference-problems.csv, with session information in
# benchmarks/results/session-info.txt.
#
# Run from the package root:  Rscript benchmarks/reference-problems.R [sizes]
# Default sizes are 8 16 32; timings depend on the machine and its load.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
sizes <- if (length(args)) as.integer(args) else c(8L, 16L, 32L)
dir.create("benchmarks/results", showWarnings = FALSE)

u_sin <- function(x, y) sin(pi * x) * sin(pi * y)
g_sin <- function(x, y) pi * cbind(cos(pi * x) * sin(pi * y), sin(pi * x) * cos(pi * y))
f_sin <- function(x, y) 2 * pi^2 * u_sin(x, y)
u_exp <- function(x, y) exp(x) * sin(y)
g_exp <- function(x, y) cbind(exp(x) * sin(y), exp(x) * cos(y))

# Each case returns list(error, ndof, assembly, solve, solver, tol) for mesh size n.
stationary <- function(make, exact, gradient) function(n) {
  sol <- solve_pde(make(n))
  d <- sol$diagnostics
  list(error = as.numeric(errornorm(sol, exact, "H1_semi", gradient = gradient)),
       error_L2 = as.numeric(errornorm(sol, exact, "L2")), ndof = d$ndof,
       assembly = d$assembly_time, solve = d$solve_time, solver = d$solver$method, tol = NA)
}

cases <- list(
  "1 Poisson, unit square" = stationary(function(n) poisson_problem(pde_mesh_rectangle(nx = n),
    source = f_sin, boundary = dirichlet_bc(0)), u_sin, g_sin),
  "2 Poisson, nonhomogeneous Dirichlet" = stationary(function(n) poisson_problem(pde_mesh_rectangle(nx = n),
    source = 0, boundary = dirichlet_bc(u_exp)), u_exp, g_exp),
  "3 Poisson, mixed conditions" = stationary(function(n) poisson_problem(pde_mesh_rectangle(nx = n),
    source = f_sin, boundary = list(dirichlet_bc(0, c("left", "right", "bottom")),
                                    neumann_bc(function(x, y) -pi * sin(pi * x), "top"))), u_sin, g_sin),
  "4 Variable diffusion" = stationary(function(n) elliptic_problem(pde_mesh_rectangle(nx = n),
    diffusion = function(x, y) 1 + x * y,
    source = function(x, y) {
      g <- g_sin(x, y); (1 + x * y) * f_sin(x, y) - (y * g[, 1] + x * g[, 2])
    }, bcs = dirichlet_bc(0)), u_sin, g_sin),
  "5 Heat equation (Crank-Nicolson, dt = h)" = function(n) {
    ex <- function(x, y, t) exp(-2 * pi^2 * t) * u_sin(x, y)
    p <- heat_problem(pde_mesh_rectangle(nx = n), initial = u_sin, bcs = dirichlet_bc(0), time = c(0, 0.1))
    sol <- solve_pde(p, dt = 0.1 / n, output_times = 0.1)
    u <- solution_at(sol, 0.1); d <- sol$diagnostics
    list(error = as.numeric(errornorm(u, function(x, y) ex(x, y, 0.1), "H1_semi",
                                      gradient = function(x, y) exp(-2 * pi^2 * 0.1) * g_sin(x, y))),
         error_L2 = as.numeric(errornorm(u, function(x, y) ex(x, y, 0.1))), ndof = d$ndof,
         assembly = d$assembly_time, solve = d$solve_time, solver = d$solver$method, tol = NA)
  },
  "6 Reaction-diffusion (linear decay, implicit Euler, dt = h^2)" = function(n) {
    ex <- function(x, y) exp(-(2 * pi^2 + 1) * 0.05) * u_sin(x, y)
    p <- reaction_diffusion_problem(pde_mesh_rectangle(nx = n), reaction = function(u, x, y, t) -u,
                                    initial = list(u = u_sin), bcs = list(u = dirichlet_bc(0)),
                                    time = c(0, 0.05))
    sol <- solve_pde(p, dt = 0.05 / n^2 * 4, method = "implicit-euler", output_times = 0.05)
    u <- solution_at(sol, 0.05); d <- sol$diagnostics
    list(error = as.numeric(errornorm(u, ex, "H1_semi", gradient = function(x, y)
           exp(-(2 * pi^2 + 1) * 0.05) * g_sin(x, y))),
         error_L2 = as.numeric(errornorm(u, ex)), ndof = d$ndof, assembly = d$assembly_time,
         solve = d$solve_time, solver = "lu (Newton)", tol = 1e-10)
  },
  "7 Advection-diffusion" = stationary(function(n) elliptic_problem(pde_mesh_rectangle(nx = n),
    advection = c(1, 0.5),
    source = function(x, y) { g <- g_sin(x, y); f_sin(x, y) + g[, 1] + 0.5 * g[, 2] },
    bcs = dirichlet_bc(0)), u_sin, g_sin),
  "8 Nonlinear diffusion k(u) = 1 + u^2" = function(n) {
    k <- function(u) 1 + u^2
    f <- function(x, y) { u <- u_sin(x, y); g <- g_sin(x, y); k(u) * f_sin(x, y) - 2 * u * rowSums(g^2) }
    sol <- solve_pde(nonlinear_problem(pde_mesh_rectangle(nx = n), diffusion = k, source = f,
                                       bcs = dirichlet_bc(0)))
    d <- sol$diagnostics
    list(error = as.numeric(errornorm(sol, u_sin, "H1_semi", gradient = g_sin)),
         error_L2 = as.numeric(errornorm(sol, u_sin)), ndof = d$ndof, assembly = NA,
         solve = d$solve_time, solver = "newton", tol = 1e-10)
  },
  "9 L-shaped polygon" = function(n) {
    L <- pde_domain("polygon", vertices = cbind(c(0, 1, 1, 0.5, 0.5, 0), c(0, 0, 0.5, 0.5, 1, 1)))
    sol <- solve_pde(poisson_problem(pde_mesh_polygon(L, h = 1 / n), source = 0, boundary = dirichlet_bc(u_exp)))
    d <- sol$diagnostics
    list(error = as.numeric(errornorm(sol, u_exp, "H1_semi", gradient = g_exp)),
         error_L2 = as.numeric(errornorm(sol, u_exp)), ndof = d$ndof, assembly = d$assembly_time,
         solve = d$solve_time, solver = d$solver$method, tol = NA)
  },
  "10 Inverse problem (diffusion coefficient)" = function(n) {
    mesh <- pde_mesh_rectangle(nx = n)
    pf <- function(theta, mesh) poisson_problem(mesh, diffusion = theta[1], source = 1, boundary = dirichlet_bc(0))
    xo <- c(0.3, 0.5, 0.7, 0.5); yo <- c(0.5, 0.5, 0.4, 0.2)
    truth <- 0.5
    yv <- 1 / truth * evaluate(solve_pde(pf(1, pde_mesh_rectangle(nx = 64)))$u, xo, yo)
    t0 <- proc.time()[["elapsed"]]
    est <- estimate_parameters(pf, theta0 = 1, observations = pde_observations(xo, yo, yv, sd = 0.01),
                               mesh = mesh, lower = 0.01)
    list(error = abs(est$estimate - truth), error_L2 = NA, ndof = est$solution$diagnostics$ndof,
         assembly = NA, solve = proc.time()[["elapsed"]] - t0, solver = "L-BFGS-B, adjoint", tol = NA)
  }
)

rows <- list()
for (nm in names(cases)) {
  prev <- NULL
  for (n in sizes) {
    invisible(gc(reset = TRUE))
    r <- cases[[nm]](n)
    g <- gc(); mem <- sum(g[, ncol(g)])  # "max used" in Mb since the reset
    h <- 1 / n
    rate <- if (is.null(prev)) NA else log(prev$error / r$error) / log(prev$h / h)
    rows[[length(rows) + 1L]] <- data.frame(problem = nm, n = n, h = h, ndof = r$ndof,
      error_H1_semi = r$error, error_L2 = r$error_L2, rate = rate,
      assembly_s = r$assembly, solve_s = r$solve, peak_mb = mem,
      solver = r$solver, tol = r$tol, stringsAsFactors = FALSE)
    prev <- list(error = r$error, h = h)
    cat(sprintf("%-60s n = %3d  ndof = %6d  error = %.3e  rate = %5.2f\n", nm, n, r$ndof, r$error,
                ifelse(is.na(rate), NaN, rate)))
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, "benchmarks/results/reference-problems.csv", row.names = FALSE)
writeLines(c(utils::capture.output(utils::sessionInfo()),
             paste("pdecompR", as.character(utils::packageVersion("pdecompR")))),
           "benchmarks/results/session-info.txt")
