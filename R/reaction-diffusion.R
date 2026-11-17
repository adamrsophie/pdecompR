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
  by_species <- function(x, what) {
    if (is.null(names(x))) {
      if (length(x) != ns) stopf("Give one %s per species.", what)
      names(x) <- species
    } else if (!setequal(names(x), species) || anyDuplicated(names(x))) {
      stopf("The names of `%s` (%s) must match the species (%s).", what,
            paste(names(x), collapse = ", "), paste(species, collapse = ", "))
    }
    x[species]
  }
  if (is.function(diffusion) || inherits(diffusion, "fe_function") || is.matrix(diffusion)) {
    diffusion <- rep(list(diffusion), ns)
  }
  if (!is.list(diffusion)) {
    if (length(diffusion) != 1L && length(diffusion) != ns) stopf("Give one diffusion coefficient or one per species.")
    diffusion <- if (is.null(names(diffusion))) as.list(rep_len(diffusion, ns)) else as.list(diffusion)
  }
  diffusion <- by_species(diffusion, "diffusion")
  for (d in diffusion) {
    check_coefficient(d, "diffusion")
    if (is.numeric(d) && d < 0) stopf("Diffusion coefficients must be nonnegative.")
  }
  if (!is.list(initial) || is.function(initial) || inherits(initial, "fe_function")) {
    initial <- list(initial)
  }
  initial <- by_species(initial, "initial")
  if (inherits(bcs, "pde_bc")) bcs <- list(bcs)
  if (length(bcs) && is.null(names(bcs))) {
    if (ns == 1L) bcs <- list(bcs) else stopf("`bcs` must be a list named by species.")
    names(bcs) <- species
  }
  if (length(bcs) && !all(names(bcs) %in% species)) {
    stopf("`bcs` names (%s) must be species names (%s).", paste(names(bcs), collapse = ", "),
          paste(species, collapse = ", "))
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
  check_tol(newton_tol, "newton_tol"); newton_maxit <- check_count(newton_maxit, "newton_maxit")
  ta <- proc.time()[["elapsed"]]
  V <- fe_space(problem$mesh, problem$degree)
  n <- V$ndof; sp <- problem$species; ns <- length(sp)
  M <- assemble_matrix(V, "mass")
  K <- lapply(problem$diffusion, function(d) {
    if (is.function(d) || inherits(d, "fe_function")) {
      dq <- coef_qp(d, element_geometry(V), name = "diffusion")
      if (any(dq < 0)) stopf("Diffusion coefficients must be nonnegative.")
    }
    assemble_matrix(V, "stiffness", d)
  })
  for (s in sp) {
    if (any(vapply(problem$bcs[[s]], function(b) b$type == "robin" && is_time_function(b$alpha),
                   logical(1)))) {
      stopf("Time-dependent Robin coefficients `alpha` are not supported; the Robin value may depend on t.")
    }
    nat <- natural_bc_terms(V, problem$bcs[[s]], problem$time[1])
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
    if (!all(is.finite(r))) stop_nonfinite("`reaction` returned non-finite values at t = %s.", fmt_num(t))
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
  newton_res <- numeric(0)
  max_rel <- 0
  bad <- NULL
  if (method == "imex-euler") {
    solvers <- list(); free <- list(); Lfd <- list()
    for (s in sp) {
      L <- M + dt * K[[s]]
      free[[s]] <- setdiff(seq_len(n), dir[[s]]$dofs)
      solvers[[s]] <- linear_solver(L[free[[s]], free[[s]], drop = FALSE], solver, ...)
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
        sol <- solvers[[s]]$solve(r, x0 = U[f, s])
        max_rel <- max(max_rel, sol$info$relative_residual)
        if (!sol$info$status %in% c("solved", "converged")) bad <- sol$info
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
        ls <- solve_linear_system(Jac, -g, method = if (solver == "auto") "lu" else solver, ...)
        if (!ls$info$status %in% c("solved", "converged")) bad <- ls$info
        max_rel <- max(max_rel, ls$info$relative_residual)
        step <- ls$x
        if (!all(is.finite(step))) break
        lam <- 1; accepted <- FALSE
        while (lam >= 1e-4) {
          Zt <- Z; zt <- as.vector(Zt); zt[freeb] <- zt[freeb] + lam * step
          Zt <- matrix(zt, n, ns, dimnames = list(NULL, sp))
          gt <- trial_eval(function() resid(Zt))
          gtn <- if (is.null(gt)) Inf else sqrt(sum(gt^2))
          if (is.finite(gtn) && gtn < (1 - 1e-4 * lam) * gn) { accepted <- TRUE; break }
          lam <- lam / 2
        }
        if (!accepted) {
          warnf("Newton line search failed at t = %s (residual %s).", fmt_num(tn), fmt_num(gn))
          break
        }
        Z <- Zt; g <- resid(Z); gn <- sqrt(sum(g^2))
      }
      if (gn > max(newton_tol * g0, 1e-12) * 10) {
        warnf("Newton did not converge at t = %s (residual %s).", fmt_num(tn), fmt_num(gn))
      }
      newton_its <- c(newton_its, it)
      newton_res <- c(newton_res, gn / max(g0, 1e-300))
      U <- Z
      if (k_out <= length(out_idx) && out_idx[k_out] == k + 1L) {
        for (s in sp) store[[s]][, k_out] <- U[, s]; k_out <- k_out + 1L
      }
    }
  }
  if (!is.null(bad)) warn_solver(bad, "linear solve of at least one time step")
  nfix <- sum(vapply(sp, function(s) length(dir[[s]]$dofs), integer(1)))
  diag <- structure(list(ndof = n * ns, free = n * ns - nfix, steps = N, dt = dt, method = method,
                         species = sp, solver = list(method = solver),
                         max_relative_residual = max_rel,
                         newton_iterations = if (length(newton_its)) newton_its else NULL,
                         newton_relative_residual = if (length(newton_res)) newton_res else NULL,
                         integral = { one <- assemble_vector(V, 1)
                                      sapply(store, function(S) as.numeric(crossprod(one, S))) },
                         assembly_time = assembly_time,
                         solve_time = proc.time()[["elapsed"]] - ts),
                    class = "pde_diagnostics")
  structure(list(times = grid[out_idx], values = store, space = V, problem = problem,
                 diagnostics = diag), class = "pde_transient")
}
