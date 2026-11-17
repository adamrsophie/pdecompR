#' Linear elasticity
#'
#' Small-strain linear elasticity in two dimensions: find the displacement
#' `u = (u_x, u_y)` with
#'
#' `-div sigma(u) = f`, `sigma(u) = lambda div(u) I + 2 mu eps(u)`,
#' `eps(u) = (grad u + grad u^T) / 2`,
#'
#' with prescribed displacements on some boundary parts and prescribed
#' tractions `sigma(u) n = g` on others. Edges without a condition are
#' traction free. The displacement uses a vector-valued continuous Lagrange
#' space: each component lies in the scalar P1 or P2 space of the mesh, and
#' the global unknowns are ordered with all `x` components first.
#'
#' The weak form is `int sigma(u) : eps(v) = int f . v + int_N g . v`; its
#' matrix is symmetric positive definite once the displacement conditions
#' remove the rigid motions (two translations and one rotation), so sparse
#' Cholesky applies. A problem without enough displacement conditions to
#' remove the rigid motions is rejected.
#'
#' Lamé parameters follow from Young's modulus `E` and Poisson's ratio `nu`:
#' `mu = E / (2 (1 + nu))` and `lambda = E nu / ((1 + nu)(1 - 2 nu))` for
#' plane strain, or `lambda = E nu / (1 - nu^2)` for plane stress. Low-order
#' elements lock for nearly incompressible materials (`nu` close to 1/2 in
#' plane strain); P2 elements reduce but do not remove the effect.
#'
#' `formulation = "mixed"` removes locking with the displacement-pressure
#' formulation: the pressure `p = lambda div u` becomes a separate unknown,
#' `sigma = 2 mu eps(u) + p I`, and the weak form is
#' `int 2 mu eps(u) : eps(v) + int p div v = int f . v + int_N g . v`,
#' `int q div u - int p q / lambda = 0`, discretized with P2 displacements
#' and P1 pressure (the Taylor-Hood pair, stable as `lambda` grows without
#' bound). Its accuracy does not degrade as `nu` approaches 1/2. The
#' solution then also has a `pressure` component. The saddle-point system is
#' solved by sparse LU or by GMRES with a block-diagonal preconditioner,
#' as for [stokes_problem()]. [solve_eigen()] and
#' [elastodynamics_problem()] need the displacement formulation.
#'
#' @param mesh A [pde_mesh()].
#' @param E,nu Young's modulus and Poisson's ratio: numbers or
#'   `function(x, y)`. `E > 0` and `-1 < nu < 1/2` (plane strain) or
#'   `-1 < nu < 1` (plane stress).
#' @param body_force Body force `f`: a length-2 vector or a `function(x, y)`
#'   returning a two-column matrix.
#' @param bcs A list of [displacement_bc()] and [traction_bc()] conditions.
#' @param degree Polynomial degree of each displacement component (1 or 2).
#' @param plane `"strain"` or `"stress"`.
#' @param quad_degree Quadrature degree (default `2 * degree + 2`).
#' @param formulation `"displacement"` or `"mixed"` (displacement-pressure,
#'   requires `degree = 2`).
#' @return An object of class `elasticity_problem`; [solve_pde()] returns a
#'   `pde_vector_solution` with the components `ux` and `uy`.
#' @examples
#' mesh <- pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 16, ny = 4)
#' prob <- elasticity_problem(mesh, E = 100, nu = 0.3,
#'                            bcs = list(displacement_bc(c(0, 0), "left"),
#'                                       traction_bc(c(0, -1), "right")))
#' sol <- solve_pde(prob)
#' sol
#' plot(sol)
#' @export
elasticity_problem <- function(mesh, E = 1, nu = 0.3, body_force = c(0, 0), bcs = list(),
                               degree = 1, plane = c("strain", "stress"), quad_degree = NULL,
                               formulation = c("displacement", "mixed")) {
  stopifnot(inherits(mesh, "pde_mesh"))
  plane <- match.arg(plane)
  formulation <- match.arg(formulation)
  if (!degree %in% 1:2) stopf("`degree` must be 1 or 2.")
  if (formulation == "mixed" && degree != 2) stopf("The mixed formulation uses Taylor-Hood elements: set `degree = 2`.")
  if (inherits(bcs, "pde_vector_bc")) bcs <- list(bcs)
  if (!all(vapply(bcs, inherits, logical(1), "pde_vector_bc"))) {
    stopf("`bcs` must be a list of displacement_bc() and traction_bc() conditions.")
  }
  check_bc_markers(bcs, mesh)
  structure(list(mesh = mesh, E = E, nu = nu, body_force = body_force, bcs = bcs,
                 degree = as.integer(degree), plane = plane, quad_degree = quad_degree,
                 formulation = formulation),
            class = c("elasticity_problem", "pde_problem"))
}

#' Boundary conditions for linear elasticity
#'
#' `displacement_bc()` prescribes displacement components on marked edges;
#' `components = "x"` fixes only `u_x` (a roller), for example.
#' `traction_bc()` prescribes the surface force `sigma(u) n`.
#'
#' @param value A length-2 vector or a `function(x, y)` returning a
#'   two-column matrix.
#' @param marker Boundary markers; `"boundary"` selects the whole boundary.
#' @param components Displacement components to prescribe (`"x"`, `"y"`).
#' @return An object of class `pde_vector_bc`.
#' @examples
#' displacement_bc(c(0, 0), "left")
#' displacement_bc(0, "bottom", components = "y")
#' traction_bc(function(x, y) cbind(0, -x), "top")
#' @export
displacement_bc <- function(value = c(0, 0), marker = "boundary", components = c("x", "y")) {
  components <- match.arg(components, c("x", "y"), several.ok = TRUE)
  if (is.numeric(value) && length(value) == 1L) value <- c(value, value)
  structure(list(type = "displacement", value = check_vector_value(value), marker = marker,
                 components = components), class = "pde_vector_bc")
}

#' @rdname displacement_bc
#' @export
traction_bc <- function(value = c(0, 0), marker = "boundary") {
  structure(list(type = "traction", value = check_vector_value(value), marker = marker),
            class = "pde_vector_bc")
}

#' @export
print.pde_vector_bc <- function(x, ...) {
  v <- if (is.function(x$value)) "function(x, y)" else paste0("(", paste(fmt_num(x$value), collapse = ", "), ")")
  comp <- if (x$type == "displacement") paste0(" [", paste(x$components, collapse = ", "), "]") else ""
  cat("<pde_vector_bc> ", x$type, comp, " = ", v, " on ", paste(x$marker, collapse = ", "), "\n", sep = "")
  invisible(x)
}

check_vector_value <- function(v) {
  if (is.function(v) || (is.numeric(v) && length(v) == 2L && all(is.finite(v)))) return(v)
  stopf("A vector boundary value must be a finite length-2 vector or a function of (x, y).")
}

# Evaluate a vector datum at points: n by 2 matrix.
eval_vector <- function(v, x, y, name) {
  if (is.numeric(v)) return(matrix(v, length(x), 2L, byrow = TRUE))
  val <- as.matrix(v(x, y))
  if (nrow(val) == 1L && length(x) > 1L) val <- matrix(val, length(x), 2L, byrow = TRUE)
  if (nrow(val) != length(x) || ncol(val) != 2L || !is.numeric(val)) {
    stopf("`%s` must return a two-column matrix with one row per point.", name)
  }
  if (any(!is.finite(val))) {
    bad <- which(!is.finite(val), arr.ind = TRUE)[1, ]
    stopf("`%s` returned a non-finite value at (%s, %s).", name, fmt_num(x[bad[1]]), fmt_num(y[bad[1]]))
  }
  val
}

lame_parameters <- function(E, nu, plane) {
  if (any(E <= 0)) stopf("Young's modulus `E` must be positive.")
  top <- if (plane == "strain") 0.5 else 1
  if (any(nu <= -1 | nu >= top)) {
    stopf("Poisson's ratio must lie in (-1, %s) for plane %s.", fmt_num(top), plane)
  }
  mu <- E / (2 * (1 + nu))
  lambda <- if (plane == "strain") E * nu / ((1 + nu) * (1 - 2 * nu)) else E * nu / (1 - nu^2)
  list(lambda = lambda, mu = mu)
}

#' @rdname discretize
#' @export
discretize.elasticity_problem <- function(problem, ...) {
  t0 <- proc.time()[["elapsed"]]
  V <- fe_space(problem$mesh, problem$degree)
  mixed <- identical(problem$formulation, "mixed")
  qd <- if (is.null(problem$quad_degree)) 2L * problem$degree + 2L else problem$quad_degree
  geo <- element_geometry(V, qd)
  n <- V$ndof
  Eq <- coef_qp(problem$E, geo, name = "E")
  nuq <- coef_qp(problem$nu, geo, name = "nu")
  lm <- lame_parameters(Eq, nuq, problem$plane)
  Z <- 0 * lm$mu
  mixed_blocks <- NULL
  if (mixed) {
    # Pressure p = lambda div u: the lambda term leaves the displacement
    # block and enters as -int p q / lambda.
    Q <- fe_space(problem$mesh, 1); gq <- element_geometry(Q, qd)
    mixed_blocks <- list(Q = Q, B = -divergence_coupling(V, Q, geo, gq),
                         C = -global_matrix(Q, local_bilinear(gq, mass = 1 / lm$lambda)),
                         Mp = global_matrix(Q, local_bilinear(gq, mass = 1 / lm$lambda + 1 / (2 * lm$mu))))
    lm$lambda <- Z
  }
  # block (b, a): test component b, trial component a,
  # D = lambda e_b e_a^T + mu delta_ab I + mu e_a e_b^T
  tens <- function(b, a) {
    D <- list(xx = Z, xy = Z, yx = Z, yy = Z)
    key <- function(p, q) paste0(c("x", "y")[p], c("x", "y")[q])
    D[[key(b, a)]] <- D[[key(b, a)]] + lm$lambda
    D[[key(a, b)]] <- D[[key(a, b)]] + lm$mu
    if (a == b) { D$xx <- D$xx + lm$mu; D$yy <- D$yy + lm$mu }
    D
  }
  # One triplet assembly for all four blocks (no block concatenation).
  D <- V$cell_dofs; nl <- ncol(D)
  I <- J <- X <- vector("list", 4L); k <- 0L
  for (bb in 1:2) for (aa in 1:2) {
    k <- k + 1L
    I[[k]] <- (bb - 1L) * n + D[, rep(seq_len(nl), times = nl)]
    J[[k]] <- (aa - 1L) * n + D[, rep(seq_len(nl), each = nl)]
    X[[k]] <- local_bilinear(geo, diff = tens(bb, aa))
  }
  A <- Matrix::sparseMatrix(i = unlist(lapply(I, as.vector)), j = unlist(lapply(J, as.vector)),
                            x = unlist(lapply(X, as.vector)), dims = c(2L * n, 2L * n))
  fq <- vec_qp(problem$body_force, geo, name = "body_force")
  b <- c(global_vector(V, local_linear(geo, f = fq$x)), global_vector(V, local_linear(geo, f = fq$y)))
  vals <- numeric(2L * n)
  fixed <- logical(2L * n)
  for (bc in problem$bcs) {
    if (bc$type == "traction") {
      comp <- function(k) function(x, y) eval_vector(bc$value, x, y, "traction")[, k]
      b <- b + c(assemble_boundary_vector(V, comp(1), bc$marker),
                 assemble_boundary_vector(V, comp(2), bc$marker))
    } else {
      d <- boundary_dofs(V, bc$marker)
      g <- eval_vector(bc$value, V$dof_coords[d, 1], V$dof_coords[d, 2], "displacement")
      if ("x" %in% bc$components) { vals[d] <- g[, 1]; fixed[d] <- TRUE }
      if ("y" %in% bc$components) { vals[n + d] <- g[, 2]; fixed[n + d] <- TRUE }
    }
  }
  dofs <- which(fixed)
  structure(list(A = A, b = b, space = V, dirichlet = list(dofs = dofs, values = vals[dofs]),
                 free = setdiff(seq_len(2L * n), dofs), assembly_time = proc.time()[["elapsed"]] - t0,
                 problem = problem, mixed = mixed_blocks), class = "pde_discretization")
}

#' @rdname solve_pde
#' @export
solve_pde.elasticity_problem <- function(problem, method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
                                         ...) {
  method <- match.arg(method)
  disc <- discretize(problem)
  V <- disc$space; n <- V$ndof
  dd <- disc$dirichlet; free <- disc$free
  # The constraints must remove the three rigid motions of every connected
  # component: the columns (1, 0), (0, 1), (-y, x) restricted to the
  # component's constrained dofs need rank 3 (coordinates centered and scaled
  # so the test does not depend on the length unit).
  xy <- rbind(V$dof_coords, V$dof_coords)
  comp_x <- rep(c(TRUE, FALSE), each = n)
  cc <- mesh_components(V$mesh)
  dof_comp <- dof_components(V, cc)
  for (k in seq_len(cc$n)) {
    on <- rep(dof_comp == k, 2L)
    ctr <- colMeans(V$dof_coords[dof_comp == k, , drop = FALSE])
    sc <- max(abs(sweep(V$dof_coords[dof_comp == k, , drop = FALSE], 2, ctr)), 1e-300)
    rx <- (xy[, 1] - ctr[1]) / sc; ry <- (xy[, 2] - ctr[2]) / sc
    Rm <- cbind(comp_x, !comp_x, ifelse(comp_x, -ry, rx))[intersect(dd$dofs, which(on)), , drop = FALSE]
    if (!nrow(Rm) || qr(Rm, tol = 1e-10)$rank < 3L) {
      stopf(paste0("The elasticity problem is singular: the displacement conditions must remove ",
                   "the rigid motions (two translations and a rotation)%s."),
            if (cc$n > 1L) sprintf(" of every connected component (component %d)", k) else "")
    }
  }
  u <- numeric(2L * n)
  u[dd$dofs] <- dd$values
  t0 <- proc.time()[["elapsed"]]
  Aff <- disc$A[free, free, drop = FALSE]
  rhs <- disc$b[free]
  if (length(dd$dofs)) rhs <- rhs - as.numeric(disc$A[free, dd$dofs, drop = FALSE] %*% dd$values)
  mx <- disc$mixed
  pressure <- NULL
  if (is.null(mx)) {
    s <- solve_linear_system(Aff, rhs, method = method, ...)
    nnz <- length(disc$A@x)
  } else {
    if (!method %in% c("auto", "lu", "gmres")) stopf("The mixed formulation is solved with method \"lu\" or \"gmres\".")
    B <- mx$B
    rhs_p <- if (length(dd$dofs)) -as.numeric(B[, dd$dofs, drop = FALSE] %*% dd$values) else numeric(nrow(B))
    s <- solve_saddle(Aff, B[, free, drop = FALSE], mx$C, rhs, rhs_p, pressure_mass = mx$Mp,
                      method = if (method == "gmres") "gmres" else "lu", ...)
    nnz <- s$nnz
    pressure <- fe_function(mx$Q, s$x[length(free) + seq_len(mx$Q$ndof)], "p")
    s$x <- s$x[seq_along(free)]
  }
  warn_solver(s$info)
  u[free] <- s$x
  diag <- structure(list(ndof = 2L * n + if (is.null(mx)) 0L else mx$Q$ndof, free = length(free), nnz = nnz,
                         assembly_time = disc$assembly_time,
                         solve_time = proc.time()[["elapsed"]] - t0, solver = s$info),
                    class = "pde_diagnostics")
  structure(list(ux = fe_function(V, u[seq_len(n)], "u_x"), uy = fe_function(V, u[n + seq_len(n)], "u_y"),
                 pressure = pressure, space = V, problem = problem, diagnostics = diag, discretization = disc),
            class = "pde_vector_solution")
}

#' @export
print.pde_vector_solution <- function(x, ...) {
  d <- x$diagnostics
  mag <- sqrt(x$ux$values^2 + x$uy$values^2)
  cat("<pde_vector_solution> P", x$space$degree, " displacement, ", d$ndof, " dofs (", d$free,
      " free); solver ", d$solver$method, " (", d$solver$status, "), relative residual ",
      fmt_num(d$solver$relative_residual, 3), "\n", sep = "")
  cat("  max |u| = ", fmt_num(max(mag)), "\n", sep = "")
  invisible(x)
}

#' @export
summary.pde_vector_solution <- function(object, ...) {
  print(object)
  invisible(object)
}

#' @rdname pde_plots
#' @param scale Displacement magnification for the deformed mesh (default:
#'   10 percent of the domain size over the largest displacement).
#' @export
plot.pde_vector_solution <- function(x, scale = NULL, main = "Deformed mesh (color: |u|)", ...) {
  P <- x$space$mesh$nodes; T <- x$space$mesh$triangles; nv <- nrow(P)
  ux <- x$ux$values[seq_len(nv)]; uy <- x$uy$values[seq_len(nv)]
  mag <- sqrt(ux^2 + uy^2)
  if (is.null(scale)) scale <- if (max(mag) > 0) 0.1 * max(diff(range(P[, 1])), diff(range(P[, 2]))) / max(mag) else 1
  Q <- P + scale * cbind(ux, uy)
  col <- grDevices::hcl.colors(64, "viridis")
  tv <- (mag[T[, 1]] + mag[T[, 2]] + mag[T[, 3]]) / 3
  k <- if (diff(range(tv)) > 0) pmin(64L, 1L + floor((tv - min(tv)) / diff(range(tv)) * 64)) else rep(1L, length(tv))
  plot(range(c(P[, 1], Q[, 1])), range(c(P[, 2], Q[, 2])), type = "n", asp = 1, xlab = "x", ylab = "y",
       main = main, ...)
  bd <- x$space$mesh$boundary
  segments(P[bd$n1, 1], P[bd$n1, 2], P[bd$n2, 1], P[bd$n2, 2], col = "grey60", lty = 2)
  polygon(as.vector(rbind(Q[T[, 1], 1], Q[T[, 2], 1], Q[T[, 3], 1], NA)),
          as.vector(rbind(Q[T[, 1], 2], Q[T[, 2], 2], Q[T[, 3], 2], NA)), col = col[k], border = "grey30",
          lwd = 0.3)
  invisible(x)
}

# Error norms of vector solutions: the exact displacement returns two columns,
# the gradient is a list of two gradient functions (one per component).
vector_errornorm <- function(u, exact, type, gradient, quad_degree = NULL, elementwise = FALSE) {
  comp <- function(k) function(x, y) eval_vector(exact, x, y, "exact")[, k]
  if (!is.null(gradient) && !(is.list(gradient) && length(gradient) == 2L)) {
    stopf("For vector solutions `gradient` must be a list of two gradient functions.")
  }
  ex <- errornorm(u$ux, comp(1), type, gradient = gradient[[1]], quad_degree = quad_degree,
                  elementwise = elementwise)
  ey <- errornorm(u$uy, comp(2), type, gradient = gradient[[2]], quad_degree = quad_degree,
                  elementwise = elementwise)
  out <- sqrt(as.numeric(ex)^2 + as.numeric(ey)^2)
  attr(out, "norm") <- type
  attr(out, "gradient") <- attr(ex, "gradient")
  if (elementwise) attr(out, "elements") <- attr(ex, "elements") + attr(ey, "elements")
  out
}

#' Stresses and strains of an elasticity solution
#'
#' Computes the strain `eps = (grad u + grad u^T) / 2`, the stress
#' `sigma = lambda tr(eps) I + 2 mu eps` and the von Mises stress
#' `sqrt(((s_xx - s_yy)^2 + (s_yy - s_zz)^2 + (s_zz - s_xx)^2) / 2 + 3 s_xy^2)`
#' of a solution of [elasticity_problem()]. The out-of-plane stress is
#' `s_zz = nu (s_xx + s_yy)` in plane strain and zero in plane stress.
#'
#' Displacement gradients are discontinuous across elements. With
#' `at = "elements"` the fields are evaluated at element centroids (exact for
#' P1, where they are constant per element). With `at = "nodes"` the centroid
#' values are averaged over the elements around each vertex with area weights
#' and returned as P1 functions, which suits plotting and VTK output; this
#' recovery is not a higher-order stress approximation.
#'
#' @param solution A `pde_vector_solution` from [solve_pde()] on an
#'   [elasticity_problem()].
#' @param at `"elements"` or `"nodes"`.
#' @return For `"elements"`, a data frame with one row per element and
#'   columns `x`, `y` (centroid), `exx`, `eyy`, `exy`, `sxx`, `syy`, `sxy`,
#'   `szz` and `von_mises`. For `"nodes"`, a named list of P1
#'   [fe_function()] objects with the same fields.
#' @examples
#' mesh <- pde_mesh_rectangle(c(0, 4), c(0, 1), nx = 16, ny = 4)
#' prob <- elasticity_problem(mesh, E = 100, nu = 0.3,
#'                            bcs = list(displacement_bc(0, "left"),
#'                                       traction_bc(c(0, -1), "right")))
#' st <- elasticity_stress(solve_pde(prob))
#' summary(st$von_mises)
#' @export
elasticity_stress <- function(solution, at = c("elements", "nodes")) {
  at <- match.arg(at)
  if (!inherits(solution, "pde_vector_solution")) stopf("`solution` must come from an elasticity problem.")
  prob <- solution$problem
  mesh <- solution$space$mesh; P <- mesh$nodes; T <- mesh$triangles
  m <- nrow(T); e <- seq_len(m)
  cx <- (P[T[, 1], 1] + P[T[, 2], 1] + P[T[, 3], 1]) / 3
  cy <- (P[T[, 1], 2] + P[T[, 2], 2] + P[T[, 3], 2]) / 3
  gx <- grad_in_elements(solution$ux, e, cx, cy)
  gy <- grad_in_elements(solution$uy, e, cx, cy)
  E <- eval_spatial(prob$E, cx, cy, name = "E")
  nu <- eval_spatial(prob$nu, cx, cy, name = "nu")
  lm <- lame_parameters(E, nu, prob$plane)
  exx <- gx$dx; eyy <- gy$dy; exy <- (gx$dy + gy$dx) / 2
  tr <- exx + eyy
  # Mixed formulation: the volumetric stress is the pressure unknown.
  vol <- if (is.null(solution$pressure)) lm$lambda * tr else grad_in_elements(solution$pressure, e, cx, cy)$value
  sxx <- vol + 2 * lm$mu * exx
  syy <- vol + 2 * lm$mu * eyy
  sxy <- 2 * lm$mu * exy
  szz <- if (prob$plane != "strain") 0 * sxx else if (is.null(solution$pressure)) nu * (sxx + syy) else vol
  vm <- sqrt(((sxx - syy)^2 + (syy - szz)^2 + (szz - sxx)^2) / 2 + 3 * sxy^2)
  out <- data.frame(x = cx, y = cy, exx = exx, eyy = eyy, exy = exy, sxx = sxx, syy = syy,
                    sxy = sxy, szz = szz, von_mises = vm)
  if (at == "elements") return(out)
  area <- abs(tri_det(P, T)) / 2
  w <- tabulate_sum(as.vector(T), rep(area, 3), nrow(P))
  V1 <- fe_space(mesh, 1)
  lapply(stats::setNames(names(out)[-(1:2)], names(out)[-(1:2)]), function(nm) {
    v <- tabulate_sum(as.vector(T), rep(out[[nm]] * area, 3), nrow(P)) / w
    fe_function(V1, v, nm)
  })
}
