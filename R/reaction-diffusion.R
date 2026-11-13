#' Reaction-diffusion systems
#'
#' `reaction_diffusion_problem()` defines a system of `S` species
#'
#' `du_s/dt - div(D_s grad u_s) = R_s(u_1, ..., u_S, x, y, t)`, `s = 1..S`,
#'
#' with per-species diffusion, initial conditions and boundary conditions.
#'
#' The reaction term is discretized by the product (group finite-element)
#' approximation: `R` is evaluated at the degrees of freedom, interpolated,
#' and integrated against the test functions, giving `M R(U)`. This is exact
#' for linear reactions and second-order accurate in space for smooth
#' nonlinear ones.
#'
#' Time integration (argument `method` of [solve_pde()]):
#'
#' * `"imex-euler"`: implicit diffusion, explicit reaction (first order);
#'   each species matrix is factored once.
#' * `"implicit-euler"`: fully implicit backward Euler with Newton's method
#'   on the coupled system and a backtracking line search. The reaction
#'   Jacobian is `jacobian` if given, otherwise computed by finite
#'   differences species by species (pointwise, so it is cheap and exact up
#'   to the difference error).
#'
#' Positivity is not guaranteed by either scheme for general reactions.
#'
#' @param mesh A [pde_mesh()].
#' @param species Character vector of species names.
#' @param diffusion Diffusion coefficients: a number, or a list or vector with
#'   one entry (number or `function(x, y)`) per species.
#' @param reaction `function(u, x, y, t)` where `u` is a matrix with one
#'   column per species (named) and one row per degree of freedom; it must
#'   return a matrix of the same shape.
#' @param jacobian Optional `function(u, x, y, t)` returning an array of
#'   dimension `c(n, S, S)` with entries `dR_a / du_b`.
#' @param initial A list with one initial condition per species.
#' @param bcs A named list of boundary conditions per species; species
#'   without an entry get zero-flux conditions.
#' @param time Start and end time.
#' @param degree Lagrange degree.
#' @return An object of class `reaction_diffusion_problem`.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 6)
#' prob <- reaction_diffusion_problem(
#'   mesh, species = "u", diffusion = 0.1,
#'   reaction = function(u, x, y, t) u * (1 - u),
#'   initial = list(u = function(x, y) 0.1 + 0.8 * x), time = c(0, 1))
#' sol <- solve_pde(prob, dt = 0.1)
#' @export
reaction_diffusion_problem <- function(mesh, species = "u", diffusion = 1, reaction,
                                       jacobian = NULL, initial, bcs = list(),
                                       time = c(0, 1), degree = 1) {
  stopifnot(inherits(mesh, "pde_mesh"))
  ns <- length(species)
  if (!is.character(species) || !ns || anyDuplicated(species)) {
    stopf("`species` must be distinct names.")
  }
  if (!is.function(reaction)) stopf("`reaction` must be a function(u, x, y, t).")
  if (!is.list(diffusion)) diffusion <- as.list(rep_len(diffusion, ns))
  if (length(diffusion) != ns) stopf("Give one diffusion coefficient per species.")
  for (d in diffusion) check_coefficient(d, "diffusion")
  names(diffusion) <- species
  if (!is.list(initial) || is.function(initial)) initial <- list(initial)
  if (length(initial) != ns) stopf("Give one initial condition per species.")
  if (is.null(names(initial))) names(initial) <- species
  initial <- initial[species]
  if (inherits(bcs, "pde_bc")) bcs <- list(bcs)
  if (length(bcs) && is.null(names(bcs))) {
    if (ns == 1L) bcs <- list(bcs) else stopf("`bcs` must be a list named by species.")
    names(bcs) <- species
  }
  bcs <- lapply(stats::setNames(species, species), function(s) normalize_bcs(bcs[[s]]))
  for (b in bcs) check_bc_markers(b, mesh)
  if (!is.numeric(time) || length(time) != 2L || time[2] <= time[1]) {
    stopf("`time` must be two increasing numbers.")
  }
  structure(list(mesh = mesh, species = species, diffusion = diffusion, reaction = reaction,
                 jacobian = jacobian, initial = initial, bcs = bcs, time = as.numeric(time),
                 degree = check_count(degree, "degree")),
            class = c("reaction_diffusion_problem", "pde_problem"))
}

#' @rdname solve_pde
#' @param newton_tol,newton_maxit Newton tolerance (on the residual norm
#'   relative to its initial value, or absolute `1e-12`) and maximum
#'   iterations per step.
#' @export
solve_pde.reaction_diffusion_problem <- function(problem, dt, method = c("imex-euler", "implicit-euler"),
                                                 output_times = NULL, solver = "auto",
                                                 newton_tol = 1e-10, newton_maxit = 25, ...) {
  method <- match.arg(method)
  check_number(dt, "dt", positive = TRUE)
  ta <- proc.time()[["elapsed"]]
  V <- fe_space(problem$mesh, problem$degree)
  n <- V$ndof; sp <- problem$species; ns <- length(sp)
  M <- assemble_matrix(V, "mass")
  K <- lapply(problem$diffusion, function(d) assemble_matrix(V, "stiffness", d))
  for (s in sp) {
    nat <- natural_bc_terms(V, problem$bcs[[s]])
    if (!is.null(nat$A)) K[[s]] <- K[[s]] + nat$A
  }
  assembly_time <- proc.time()[["elapsed"]] - ta
  t0 <- problem$time[1]; t1 <- problem$time[2]
  N <- max(1L, as.integer(ceiling((t1 - t0) / dt - 1e-10)))
  dt <- (t1 - t0) / N
  grid <- t0 + (0:N) * dt
  out_idx <- output_index(output_times, grid, dt)
  xy <- V$dof_coords
  U <- vapply(sp, function(s) initial_values(problem$initial[[s]], V), numeric(n))
  U <- matrix(U, n, ns, dimnames = list(NULL, sp))
  dir <- lapply(sp, function(s) dirichlet_data(V, problem$bcs[[s]], t0))
  names(dir) <- sp
  for (s in sp) if (length(dir[[s]]$dofs)) U[dir[[s]]$dofs, s] <- dir[[s]]$values
  reac <- function(U, t) {
    r <- as.matrix(problem$reaction(U, xy[, 1], xy[, 2], t))
    if (!all(dim(r) == dim(U))) stopf("`reaction` must return a matrix with %d rows and %d columns.", n, ns)
    dimnames(r) <- list(NULL, sp)
    r
  }
  natload <- function(s, t) natural_bc_terms(V, problem$bcs[[s]], t)$b
  store <- lapply(sp, function(s) matrix(NA_real_, n, length(out_idx)))
  names(store) <- sp
  k_out <- 1L
  if (out_idx[1] == 1L) { for (s in sp) store[[s]][, 1] <- U[, s]; k_out <- 2L }
  ts <- proc.time()[["elapsed"]]
  newton_its <- integer(0)
  max_rel <- 0
  if (method == "imex-euler") {
    solvers <- list(); free <- list(); Lfd <- list()
    for (s in sp) {
      L <- M + dt * K[[s]]
      free[[s]] <- setdiff(seq_len(n), dir[[s]]$dofs)
      solvers[[s]] <- linear_solver(L[free[[s]], free[[s]], drop = FALSE], solver)
      Lfd[[s]] <- L[free[[s]], dir[[s]]$dofs, drop = FALSE]
    }
    for (k in seq_len(N)) {
      tn <- grid[k + 1L]
      R <- reac(U, grid[k])
      Unew <- U
      for (s in sp) {
        rhs <- as.numeric(M %*% (U[, s] + dt * R[, s])) + dt * natload(s, tn)
        dd <- dirichlet_data(V, problem$bcs[[s]], tn)
        f <- free[[s]]
        r <- rhs[f]
        if (length(dd$dofs)) r <- r - as.numeric(Lfd[[s]] %*% dd$values)
        sol <- solvers[[s]]$solve(r)
        max_rel <- max(max_rel, sol$info$relative_residual)
        Unew[f, s] <- sol$x
        if (length(dd$dofs)) Unew[dd$dofs, s] <- dd$values
      }
      U <- Unew
      if (k_out <= length(out_idx) && out_idx[k_out] == k + 1L) {
        for (s in sp) store[[s]][, k_out] <- U[, s]; k_out <- k_out + 1L
      }
    }
  } else {
    Mb <- Matrix::bdiag(rep(list(M), ns))
    Kb <- Matrix::bdiag(K)
    Lb <- Mb + dt * Kb
    fixed <- unlist(lapply(seq_len(ns), function(i) (i - 1L) * n + dir[[sp[i]]]$dofs))
    freeb <- setdiff(seq_len(n * ns), fixed)
    jac_r <- function(U, t) {
      if (!is.null(problem$jacobian)) return(problem$jacobian(U, xy[, 1], xy[, 2], t))
      J <- array(0, c(n, ns, ns))
      R0 <- reac(U, t)
      for (b in seq_len(ns)) {
        h <- 1e-7 * pmax(1, abs(U[, b]))
        Up <- U; Up[, b] <- Up[, b] + h
        J[, , b] <- (reac(Up, t) - R0) / h
      }
      J
    }
    for (k in seq_len(N)) {
      tn <- grid[k + 1L]
      zold <- as.vector(U)
      Fn <- unlist(lapply(sp, function(s) natload(s, tn)))
      Z <- U
      for (s in sp) {
        dd <- dirichlet_data(V, problem$bcs[[s]], tn)
        if (length(dd$dofs)) Z[dd$dofs, s] <- dd$values
      }
      resid <- function(Z) {
        z <- as.vector(Z)
        g <- as.numeric(Lb %*% z) - as.numeric(Mb %*% zold) -
          dt * as.numeric(Mb %*% as.vector(reac(Z, tn))) - dt * Fn
        g[freeb]
      }
      g <- resid(Z); g0 <- sqrt(sum(g^2)); gn <- g0
      it <- 0L
      while (gn > max(newton_tol * g0, 1e-12) && it < newton_maxit) {
        it <- it + 1L
        J <- jac_r(Z, tn)
        blocks <- lapply(seq_len(ns), function(a) lapply(seq_len(ns), function(b) {
          Matrix::Diagonal(x = J[, a, b])
        }))
        JR <- do.call(rbind, lapply(blocks, function(row) do.call(cbind, row)))
        Jac <- (Lb - dt * Mb %*% JR)[freeb, freeb, drop = FALSE]
        step <- solve_linear_system(Jac, -g, method = "lu")$x
        lam <- 1
        repeat {
          Zt <- Z; zt <- as.vector(Zt); zt[freeb] <- zt[freeb] + lam * step
          Zt <- matrix(zt, n, ns, dimnames = list(NULL, sp))
          gt <- resid(Zt)
          if (sqrt(sum(gt^2)) < (1 - 1e-4 * lam) * gn || lam < 1e-4) break
          lam <- lam / 2
        }
        Z <- Zt; g <- gt; gn <- sqrt(sum(g^2))
      }
      if (gn > max(newton_tol * g0, 1e-12) * 10) {
        warnf("Newton did not converge at t = %s (residual %s).", fmt_num(tn), fmt_num(gn))
      }
      newton_its <- c(newton_its, it)
      U <- Z
      if (k_out <= length(out_idx) && out_idx[k_out] == k + 1L) {
        for (s in sp) store[[s]][, k_out] <- U[, s]; k_out <- k_out + 1L
      }
    }
  }
  diag <- structure(list(ndof = n, free = n, steps = N, dt = dt, method = method,
                         species = sp, solver = list(method = solver),
                         max_relative_residual = max_rel,
                         newton_iterations = if (length(newton_its)) newton_its else NULL,
                         assembly_time = assembly_time,
                         solve_time = proc.time()[["elapsed"]] - ts),
                    class = "pde_diagnostics")
  structure(list(times = grid[out_idx], values = store, space = V, problem = problem,
                 diagnostics = diag), class = "pde_transient")
}
