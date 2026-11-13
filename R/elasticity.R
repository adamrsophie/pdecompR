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
                               degree = 1, plane = c("strain", "stress"), quad_degree = NULL) {
  stopifnot(inherits(mesh, "pde_mesh"))
  plane <- match.arg(plane)
  if (!degree %in% 1:2) stopf("`degree` must be 1 or 2.")
  if (inherits(bcs, "pde_vector_bc")) bcs <- list(bcs)
  if (!all(vapply(bcs, inherits, logical(1), "pde_vector_bc"))) {
    stopf("`bcs` must be a list of displacement_bc() and traction_bc() conditions.")
  }
  check_bc_markers(bcs, mesh)
  structure(list(mesh = mesh, E = E, nu = nu, body_force = body_force, bcs = bcs,
                 degree = as.integer(degree), plane = plane, quad_degree = quad_degree),
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
  geo <- element_geometry(V, problem$quad_degree)
  n <- V$ndof
  Eq <- coef_qp(problem$E, geo, name = "E")
  nuq <- coef_qp(problem$nu, geo, name = "nu")
  lm <- lame_parameters(Eq, nuq, problem$plane)
  Z <- 0 * lm$mu
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
  blk <- function(b, a) global_matrix(V, local_bilinear(geo, diff = tens(b, a)))
  A <- rbind(cbind(blk(1, 1), blk(1, 2)), cbind(blk(2, 1), blk(2, 2)))
  fq <- vec_qp(problem$body_force, geo, name = "body_force")
  b <- c(global_vector(V, local_linear(geo, f = fq$x)), global_vector(V, local_linear(geo, f = fq$y)))
  vals <- rep(NA_real_, 2L * n)
  for (bc in problem$bcs) {
    if (bc$type == "traction") {
      comp <- function(k) function(x, y) eval_vector(bc$value, x, y, "traction")[, k]
      b <- b + c(assemble_boundary_vector(V, comp(1), bc$marker),
                 assemble_boundary_vector(V, comp(2), bc$marker))
    } else {
      d <- boundary_dofs(V, bc$marker)
      g <- eval_vector(bc$value, V$dof_coords[d, 1], V$dof_coords[d, 2], "displacement")
      if ("x" %in% bc$components) vals[d] <- g[, 1]
      if ("y" %in% bc$components) vals[n + d] <- g[, 2]
    }
  }
  dofs <- which(!is.na(vals))
  structure(list(A = A, b = b, space = V, dirichlet = list(dofs = dofs, values = vals[dofs]),
                 free = setdiff(seq_len(2L * n), dofs), assembly_time = proc.time()[["elapsed"]] - t0,
                 problem = problem), class = "pde_discretization")
}

#' @rdname solve_pde
#' @export
solve_pde.elasticity_problem <- function(problem, method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
                                         ...) {
  method <- match.arg(method)
  disc <- discretize(problem)
  V <- disc$space; n <- V$ndof
  dd <- disc$dirichlet; free <- disc$free
  # The constraints must remove the three rigid motions: the 3 columns
  # (1, 0), (0, 1), (-y, x) restricted to the constrained dofs need rank 3.
  xy <- rbind(V$dof_coords, V$dof_coords)
  comp_x <- rep(c(TRUE, FALSE), each = n)
  Rm <- cbind(comp_x, !comp_x, ifelse(comp_x, -xy[, 2], xy[, 1]))[dd$dofs, , drop = FALSE]
  if (!length(dd$dofs) || qr(Rm, tol = 1e-10)$rank < 3L) {
    stopf(paste0("The elasticity problem is singular: the displacement conditions must remove ",
                 "the rigid motions (two translations and a rotation)."))
  }
  u <- numeric(2L * n)
  u[dd$dofs] <- dd$values
  t0 <- proc.time()[["elapsed"]]
  Aff <- disc$A[free, free, drop = FALSE]
  rhs <- disc$b[free]
  if (length(dd$dofs)) rhs <- rhs - as.numeric(disc$A[free, dd$dofs, drop = FALSE] %*% dd$values)
  s <- solve_linear_system(Aff, rhs, method = method, ...)
  warn_solver(s$info)
  u[free] <- s$x
  diag <- structure(list(ndof = 2L * n, free = length(free), nnz = length(disc$A@x),
                         assembly_time = disc$assembly_time,
                         solve_time = proc.time()[["elapsed"]] - t0, solver = s$info),
                    class = "pde_diagnostics")
  structure(list(ux = fe_function(V, u[seq_len(n)], "u_x"), uy = fe_function(V, u[n + seq_len(n)], "u_y"),
                 space = V, problem = problem, diagnostics = diag, discretization = disc),
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
vector_errornorm <- function(u, exact, type, gradient) {
  comp <- function(k) function(x, y) eval_vector(exact, x, y, "exact")[, k]
  if (!is.null(gradient) && !(is.list(gradient) && length(gradient) == 2L)) {
    stopf("For vector solutions `gradient` must be a list of two gradient functions.")
  }
  ex <- errornorm(u$ux, comp(1), type, gradient = gradient[[1]])
  ey <- errornorm(u$uy, comp(2), type, gradient = gradient[[2]])
  out <- sqrt(as.numeric(ex)^2 + as.numeric(ey)^2)
  attr(out, "norm") <- type
  attr(out, "gradient") <- attr(ex, "gradient")
  out
}
