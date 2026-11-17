#' Stokes flow with Taylor-Hood elements
#'
#' `stokes_problem()` defines the incompressible Stokes equations
#'
#' `-div(nu grad u) + grad p = f`, `div u = 0`
#'
#' for velocity `u` and pressure `p`, discretized with the inf-sup stable
#' Taylor-Hood pair: continuous P2 velocity and P1 pressure. Velocity
#' conditions are given with [velocity_bc()] (any components, data may vary
#' in space); [traction_bc()] prescribes the pseudo-traction
#' `(nu grad u - p I) n`, and boundary parts with no condition get the
#' do-nothing outflow condition (zero pseudo-traction). When the velocity is
#' prescribed on the whole boundary, the pressure is fixed by a zero-mean
#' gauge.
#'
#' [solve_pde()] solves the saddle-point system by sparse LU
#' (`method = "lu"`, the default) or by right-preconditioned GMRES with the
#' block-diagonal preconditioner `diag(A, M_p / nu)` (`method = "gmres"`),
#' where `A` is the velocity Laplacian (factorized once) and `M_p` the
#' pressure mass matrix; its iteration count does not grow with mesh
#' refinement.
#'
#' @param mesh A [pde_mesh()].
#' @param viscosity Kinematic viscosity `nu > 0` (number or `function(x, y)`).
#' @param body_force Body force `f` (length-2 vector or `function(x, y)`
#'   returning a two-column matrix).
#' @param bcs List of [velocity_bc()] and [traction_bc()] conditions.
#' @return An object of class `stokes_problem`. Its solution has components
#'   `velocity` (a vector solution with `ux` and `uy`, usable with
#'   [errornorm()]), `pressure` (an [fe_function()]) and `diagnostics`
#'   (including the L2 norm of `div u_h`).
#' @references Taylor, C. and Hood, P. (1973). A numerical solution of the
#'   Navier-Stokes equations using the finite element technique. *Computers
#'   and Fluids*, 1, 73-100. Elman, H., Silvester, D. and Wathen, A. (2014).
#'   *Finite Elements and Fast Iterative Solvers*, 2nd edition. Oxford
#'   University Press.
#' @examples
#' # Lid-driven cavity
#' m <- pde_mesh_rectangle(nx = 8)
#' p <- stokes_problem(m, bcs = list(velocity_bc(c(0, 0)), velocity_bc(c(1, 0), "top")))
#' s <- solve_pde(p)
#' range(s$velocity$ux$values)
#' @export
stokes_problem <- function(mesh, viscosity = 1, body_force = c(0, 0), bcs = list()) {
  stopifnot(inherits(mesh, "pde_mesh"))
  check_coefficient(viscosity, "viscosity")
  if (is.numeric(viscosity) && viscosity <= 0) stopf("`viscosity` must be positive.")
  if (inherits(bcs, "pde_vector_bc")) bcs <- list(bcs)
  if (!all(vapply(bcs, inherits, logical(1), "pde_vector_bc"))) {
    stopf("`bcs` must be a list of velocity_bc() and traction_bc() conditions.")
  }
  check_bc_markers(bcs, mesh)
  structure(list(mesh = mesh, viscosity = viscosity, body_force = check_vector_value(body_force), bcs = bcs,
                 degree = 2L), class = c("stokes_problem", "pde_problem"))
}

#' @rdname stokes_problem
#' @param value Velocity `c(ux, uy)` or `function(x, y)` returning a two-column
#'   matrix.
#' @param marker Boundary markers.
#' @param components Velocity components to prescribe.
#' @export
velocity_bc <- function(value = c(0, 0), marker = "boundary", components = c("x", "y")) {
  displacement_bc(value, marker, components)
}

#' @rdname solve_pde
#' @export
solve_pde.stokes_problem <- function(problem, method = c("lu", "gmres"), tol = 1e-10, ...) {
  method <- match.arg(method)
  t0 <- proc.time()[["elapsed"]]
  mesh <- problem$mesh
  V <- fe_space(mesh, 2); Q <- fe_space(mesh, 1)
  n <- V$ndof; nq <- Q$ndof
  gv <- element_geometry(V, 4L); gq <- element_geometry(Q, 4L)
  nuq <- coef_qp(problem$viscosity, gv, name = "viscosity")
  if (any(nuq <= 0)) stopf("`viscosity` must be positive.")
  K <- global_matrix(V, local_bilinear(gv, diff = nuq))
  # B_ij = -int psi_i d(phi_j)/dx_c for pressure test psi_i and velocity phi_j.
  bmat <- function(c) {
    Bl <- array(0, c(gv$m, gq$nloc, gv$nloc))
    for (k in seq_len(gv$nq)) {
      wd <- gv$w[k] * gv$det
      G <- if (c == 1L) gv$gx[[k]] else gv$gy[[k]]
      for (j in seq_len(gv$nloc)) for (i in seq_len(gq$nloc)) Bl[, i, j] <- Bl[, i, j] - wd * gq$phi[k, i] * G[, j]
    }
    Matrix::sparseMatrix(i = as.vector(Q$cell_dofs[, rep(seq_len(gq$nloc), times = gv$nloc)]),
                         j = as.vector(V$cell_dofs[, rep(seq_len(gv$nloc), each = gq$nloc)]),
                         x = as.vector(Bl), dims = c(nq, n))
  }
  Bx <- bmat(1L); By <- bmat(2L)
  Z <- function(r, c) Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(r, c))
  B <- cbind(Bx, By)
  Avel <- rbind(cbind(K, Z(n, n)), cbind(Z(n, n), K))
  fq <- vec_qp(problem$body_force, gv, name = "body_force")
  b <- c(global_vector(V, local_linear(gv, f = fq$x)), global_vector(V, local_linear(gv, f = fq$y)))
  vals <- numeric(2L * n); fixed <- logical(2L * n)
  for (bc in problem$bcs) {
    if (bc$type == "traction") {
      comp <- function(k) function(x, y) eval_vector(bc$value, x, y, "traction")[, k]
      b <- b + c(assemble_boundary_vector(V, comp(1), bc$marker), assemble_boundary_vector(V, comp(2), bc$marker))
    } else {
      d <- boundary_dofs(V, bc$marker)
      g <- eval_vector(bc$value, V$dof_coords[d, 1], V$dof_coords[d, 2], "velocity")
      if ("x" %in% bc$components) { vals[d] <- g[, 1]; fixed[d] <- TRUE }
      if ("y" %in% bc$components) { vals[n + d] <- g[, 2]; fixed[n + d] <- TRUE }
    }
  }
  bd <- boundary_dofs(V)
  enclosed <- all(fixed[bd]) && all(fixed[n + bd])
  dofs <- which(fixed); fv <- which(!fixed)
  rhs_u <- b[fv] - as.numeric(Avel[fv, dofs, drop = FALSE] %*% vals[dofs])
  rhs_p <- -as.numeric(B[, dofs, drop = FALSE] %*% vals[dofs])
  Aff <- Avel[fv, fv, drop = FALSE]; Bf <- B[, fv, drop = FALSE]
  mp <- assemble_vector(Q, 1)
  ng <- if (enclosed) 1L else 0L
  S <- rbind(cbind(Aff, Matrix::t(Bf), Z(length(fv), ng)),
             cbind(Bf, Z(nq, nq), if (ng) Matrix::Matrix(mp, ncol = 1, sparse = TRUE) else Z(nq, 0)),
             cbind(Z(ng, length(fv)), if (ng) Matrix::Matrix(mp, nrow = 1, sparse = TRUE) else Z(0, nq), Z(ng, ng)))
  rhs <- c(rhs_u, rhs_p, numeric(ng))
  assembly_time <- proc.time()[["elapsed"]] - t0
  ts <- proc.time()[["elapsed"]]
  if (method == "lu") {
    s <- solve_linear_system(S, rhs, method = "lu")
  } else {
    sa <- linear_solver(Aff)
    sm <- linear_solver(assemble_matrix(Q, "mass", coef = function(x, y) 1 / eval_spatial(problem$viscosity, x, y)))
    nu_ <- length(fv)
    pc <- function(r) {
      c(sa$solve(r[seq_len(nu_)])$x, sm$solve(r[nu_ + seq_len(nq)])$x, r[nu_ + nq + seq_len(ng)])
    }
    s <- solve_linear_system(S, rhs, method = "gmres", preconditioner = pc, tol = tol, ...)
  }
  warn_solver(s$info)
  u <- vals; u[fv] <- s$x[seq_len(length(fv))]
  p <- s$x[length(fv) + seq_len(nq)]
  ux <- fe_function(V, u[seq_len(n)], "u_x"); uy <- fe_function(V, u[n + seq_len(n)], "u_y")
  # L2 norm of the discrete divergence.
  qx <- fe_qp(ux, gv); qy <- fe_qp(uy, gv)
  divn <- sqrt(sum(((qx$dx + qy$dy)^2 %*% gv$w) * gv$det))
  diag <- structure(list(ndof = 2L * n + nq, free = length(fv) + nq, nnz = length(S@x),
                         assembly_time = assembly_time, solve_time = proc.time()[["elapsed"]] - ts,
                         solver = s$info, pressure_gauge = if (enclosed) "zero mean" else "none",
                         divergence_norm = divn),
                    class = "pde_diagnostics")
  vel <- structure(list(ux = ux, uy = uy, space = V, problem = problem, diagnostics = diag),
                   class = "pde_vector_solution")
  structure(list(velocity = vel, pressure = fe_function(Q, p, "p"), problem = problem, diagnostics = diag),
            class = "pde_stokes_solution")
}

#' @export
print.pde_stokes_solution <- function(x, ...) {
  d <- x$diagnostics
  cat("<pde_stokes_solution> Taylor-Hood P2-P1, ", d$ndof, " dofs; solver ", d$solver$method, " (",
      d$solver$status, ")", if (!is.na(d$solver$iterations)) paste0(", ", d$solver$iterations, " iterations"),
      "\n", sep = "")
  cat("  ||div u_h|| = ", fmt_num(d$divergence_norm, 3), "; pressure gauge: ", d$pressure_gauge, "\n", sep = "")
  invisible(x)
}
