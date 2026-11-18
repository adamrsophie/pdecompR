#' Errors against a reference solution
#'
#' `errornorm()` computes the error between a finite-element function and an
#' exact solution by element-wise quadrature of a rule of higher degree than
#' the one used for assembly (default `2 p + 4` for degree `p`). Nodal errors
#' are not used: they measure a different (superconvergent) quantity.
#'
#' Norms:
#'
#' * `"L2"`: `||u - u_h||_{L2}`;
#' * `"H1_semi"`: `|u - u_h|_{H1} = ||grad(u - u_h)||_{L2}`;
#' * `"H1"`: the full norm `sqrt(L2^2 + H1_semi^2)`;
#' * `"energy"`: `sqrt(int k |grad(u - u_h)|^2)` with diffusion `k`.
#'
#' Gradient-based norms need the exact gradient. When `gradient` is not given
#' it is approximated by central differences with step `1e-5` times the
#' domain diameter; the result is then labeled with attribute
#' `gradient = "finite difference"`. Its accuracy depends on the smoothness of
#' `exact` (the truncation error is of order step squared times the third
#' derivative), and it is wrong near jumps of the gradient, so pass `gradient`
#' for verification work. The step does not change when the domain is
#' translated; on a domain far from the origin it is raised where the
#' coordinates cannot represent it, and the rounding error then grows like
#' `1e-16 |x|` divided by the step.
#'
#' @param u An [fe_function()], `pde_solution` or `pde_vector_solution` (from
#'   [elasticity_problem()]; then `exact` returns two columns and `gradient` is a
#'   list of two gradient functions, one per component).
#' @param exact Exact solution `function(x, y)`, or an [fe_function()]
#'   reference solution on any mesh of the same domain, possibly of another
#'   degree; its gradient is then used for gradient norms. The integral runs
#'   on a partition on which both functions are polynomials: the finer of two
#'   nested meshes (for example from [refine_uniform()] or [refine_mesh()]),
#'   or the overlay of all pairwise triangle intersections for meshes that
#'   are not nested (slower). The result is exact up to quadrature and
#'   symmetric in the two functions.
#' @param type Norm type.
#' @param gradient Exact gradient `function(x, y)` returning a two-column
#'   matrix.
#' @param diffusion Diffusion coefficient for the energy norm: a number, a
#'   function of `(x, y)` or a constant 2 by 2 matrix.
#' @param quad_degree Quadrature exactness degree.
#' @param elementwise Return element contributions (squared) as attribute
#'   `elements`.
#' @return The error (a number) with attributes.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
#' sol <- solve_pde(poisson_problem(mesh, source = f, boundary = dirichlet_bc(0)))
#' errornorm(sol, function(x, y) sin(pi * x) * sin(pi * y), "L2")
#' @export
errornorm <- function(u, exact, type = c("L2", "H1_semi", "H1", "energy"), gradient = NULL,
                      diffusion = 1, quad_degree = NULL, elementwise = FALSE) {
  type <- match.arg(type)
  if (inherits(u, "pde_vector_solution")) {
    if (type == "energy") stopf("Use \"L2\", \"H1_semi\" or \"H1\" for vector solutions.")
    return(vector_errornorm(u, exact, type, gradient, quad_degree, elementwise))
  }
  if (inherits(u, "pde_transient")) {
    stopf("Use `errornorm(solution_at(sol, time), exact)` for time-dependent solutions.")
  }
  if (inherits(u, "pde_solution")) u <- u$u
  if (!inherits(u, "fe_function")) stopf("`u` must be a finite-element function or solution.")
  label <- "exact"
  ref_degree <- 0L
  if (inherits(exact, "fe_function")) {
    # Reference solution, for example on a finer mesh. Quadrature must run on
    # a partition on which both functions are polynomials: the mesh of `u`
    # when its triangles lie inside those of the reference mesh, otherwise the
    # reference mesh (when nested the other way), with contributions summed
    # back to the triangles of `u`.
    ref <- exact
    fast <- same_space(u$space, ref$space) && is.null(gradient)
    if (fast) {
      # Same space: the error is the finite-element function u - ref, whose
      # values and gradients are read at the quadrature points directly,
      # without locating points in the reference mesh.
      u$values <- u$values - ref$values
      exact <- function(x, y) numeric(length(x))
      gradient <- function(x, y) matrix(0, length(x), 2)
      label <- "finite-element reference"
      ref_degree <- ref$space$degree
    } else if (is.null(nested_parent(u$space$mesh, ref$space$mesh))) {
      par <- nested_parent(ref$space$mesh, u$space$mesh)
      if (is.null(quad_degree)) quad_degree <- 2L * max(u$space$degree, ref$space$degree) + 4L
      if (is.null(par)) {
        # Meshes that are not nested: integrate on their overlay.
        e2 <- overlay_error(u, ref, type, diffusion, quad_degree)
      } else {
        r <- errornorm(ref, u, type, diffusion = diffusion, quad_degree = quad_degree, elementwise = TRUE)
        e2 <- numeric(nrow(u$space$mesh$triangles))
        s2 <- rowsum(attr(r, "elements"), par, reorder = TRUE)
        e2[as.integer(rownames(s2))] <- s2[, 1]
      }
      out <- sqrt(sum(e2))
      attr(out, "norm") <- type
      attr(out, "gradient") <- "finite-element reference"
      if (elementwise) attr(out, "elements") <- e2
      return(out)
    }
    if (!fast) {
      exact <- function(x, y) evaluate(ref, x, y)
      if (is.null(gradient)) {
        gradient <- function(x, y) as.matrix(evaluate(ref, x, y, gradient = TRUE)[, 2:3])
        label <- "finite-element reference"
      }
      ref_degree <- ref$space$degree
    }
  }
  V <- u$space
  if (is.null(quad_degree)) quad_degree <- 2L * max(V$degree, ref_degree) + 4L
  geo <- element_geometry(V, quad_degree)
  uh <- fe_qp(u, geo)
  xx <- as.vector(geo$x); yy <- as.vector(geo$y)
  w <- function(vals) as.vector((vals %*% geo$w) * geo$det)
  e2 <- numeric(geo$m)
  if (type %in% c("L2", "H1")) {
    ue <- matrix(eval_spatial(exact, xx, yy, name = "exact"), geo$m)
    e2 <- e2 + w((ue - uh$value)^2)
  }
  if (type %in% c("H1_semi", "H1", "energy")) {
    if (is.null(gradient)) {
      label <- "finite difference"
      bb <- apply(V$mesh$nodes, 2, range)
      gradient <- fd_gradient(exact, scale = sqrt(sum((bb[2, ] - bb[1, ])^2)))
    }
    g <- as.matrix(gradient(xx, yy))
    if (nrow(g) != length(xx) || ncol(g) != 2L) {
      stopf("`gradient` must return a two-column matrix with one row per point.")
    }
    ex <- matrix(g[, 1], geo$m) - uh$dx; ey <- matrix(g[, 2], geo$m) - uh$dy
    d2 <- ex^2 + ey^2
    if (type == "energy") {
      d2 <- if (inherits(diffusion, "pde_tensor") || (is.matrix(diffusion) && is.numeric(diffusion))) {
        energy_density(diffusion, ex, ey, geo$x, geo$y)
      } else {
        d2 * coef_qp(diffusion, geo, name = "diffusion")
      }
    }
    e2 <- e2 + w(d2)
  }
  out <- sqrt(sum(e2))
  attr(out, "norm") <- type
  attr(out, "gradient") <- label
  if (elementwise) attr(out, "elements") <- e2
  out
}

# Central-difference gradient with step h * scale, where scale is the
# diameter of the domain: the step follows the size of the domain, not the
# distance from the origin, so translating the domain does not change it.
# The step is raised where coordinates are so large that it is not
# representable, and the divisor is the step actually taken in floating
# point.
fd_gradient <- function(f, h = 1e-5, scale = 1) {
  function(x, y) {
    st <- function(v) pmax(h * scale, 4 * .Machine$double.eps * abs(v))
    hx <- st(x); hy <- st(y)
    xp <- x + hx; xm <- x - hx; yp <- y + hy; ym <- y - hy
    cbind((f(xp, y) - f(xm, y)) / (xp - xm), (f(x, yp) - f(x, ym)) / (yp - ym))
  }
}

#' Mesh convergence studies
#'
#' Solves a problem on a sequence of meshes and reports errors and observed
#' convergence orders `log(e_{i-1} / e_i) / log(h_{i-1} / h_i)`.
#'
#' When the problem's mesh came from [pde_mesh_rectangle()] and `mesh_sizes`
#' is given, meshes with `nx = mesh_sizes[i]` cells (and the same aspect) are
#' generated. Otherwise the problem's mesh is refined uniformly
#' `refinements` times.
#'
#' Observed orders on coarse meshes can deviate from the asymptotic ones
#' (pre-asymptotic behavior). For P1 and a smooth solution the asymptotic
#' orders are 2 in L2 and 1 in H1; for P2 they are 3 and 2.
#'
#' @param problem A stationary `pde_problem` (elliptic, nonlinear or
#'   elasticity). Variational problems are tied to their function space and
#'   are not supported.
#' @param exact Exact solution `function(x, y)`.
#' @param gradient Exact gradient (optional, see [errornorm()]).
#' @param mesh_sizes Cell counts for rectangle meshes.
#' @param refinements Number of uniform refinements otherwise.
#' @param norms Norms to report.
#' @param ... Passed to [solve_pde()].
#' @return A data frame of class `pde_convergence` with columns `h`, `ndof`,
#'   one column per norm and the observed rates. `summary()` reports the
#'   finest errors and the last and mean observed rates per norm.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 4)
#' f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
#' prob <- poisson_problem(mesh, source = f, boundary = dirichlet_bc(0))
#' convergence_study(prob, function(x, y) sin(pi * x) * sin(pi * y),
#'                   mesh_sizes = c(4, 8, 16))
#' @export
convergence_study <- function(problem, exact, gradient = NULL, mesh_sizes = NULL,
                              refinements = 3, norms = c("L2", "H1_semi"), ...) {
  stopifnot(inherits(problem, "pde_problem"))
  if (inherits(problem, "variational_problem")) {
    stopf("Convergence studies need a problem that can be posed on new meshes; variational problems cannot.")
  }
  if (!inherits(problem, c("elliptic_problem", "nonlinear_problem", "elasticity_problem"))) {
    stopf("Convergence studies support stationary elliptic, nonlinear and elasticity problems.")
  }
  if ("energy" %in% norms && inherits(problem, "nonlinear_problem")) {
    stopf("The energy norm is not defined for nonlinear problems; use \"L2\" or \"H1_semi\".")
  }
  args <- problem$mesh$info$args
  if (!is.null(mesh_sizes)) {
    if (isTRUE(problem$mesh$info$markers_edited)) {
      stopf("The mesh has relabeled boundary markers that new meshes would lose; use `refinements`.")
    }
    if (!identical(problem$mesh$info$generator, "rectangle")) {
      stopf("`mesh_sizes` needs a problem on a mesh from `pde_mesh_rectangle()`; use `refinements`.")
    }
    ratio <- args$ny / args$nx
    meshes <- lapply(mesh_sizes, function(n) {
      pde_mesh_rectangle(args$x, args$y, nx = n, ny = max(1, round(n * ratio)), pattern = args$pattern)
    })
  } else {
    refinements <- check_count(refinements, "refinements", min = 1L)
    meshes <- vector("list", refinements + 1L)
    meshes[[1]] <- problem$mesh
    for (i in seq_len(refinements)) meshes[[i + 1L]] <- refine_uniform(meshes[[i]])
  }
  rows <- lapply(meshes, function(mesh) {
    sol <- solve_pde(update_mesh(problem, mesh), ...)
    errs <- vapply(norms, function(nm) as.numeric(errornorm(sol, exact, nm, gradient = gradient,
                                                             diffusion = problem$diffusion)),
                   numeric(1))
    c(h = mesh_size(mesh), ndof = sol$diagnostics$ndof, errs)
  })
  out <- as.data.frame(do.call(rbind, rows))
  for (nm in norms) {
    r <- c(NA, diff(log(out[[nm]])) / diff(log(out$h)))
    out[[paste0("rate_", nm)]] <- r
  }
  structure(out, class = c("pde_convergence", "data.frame"), norms = norms)
}

#' @export
summary.pde_convergence <- function(object, ...) {
  norms <- attr(object, "norms")
  last <- nrow(object)
  out <- data.frame(norm = norms,
                    finest_error = vapply(norms, function(nm) object[[nm]][last], numeric(1)),
                    last_rate = vapply(norms, function(nm) object[[paste0("rate_", nm)]][last], numeric(1)),
                    mean_rate = vapply(norms, function(nm) mean(object[[paste0("rate_", nm)]][-1]), numeric(1)),
                    row.names = NULL)
  structure(list(rates = out, h = range(object$h), ndof = range(object$ndof)),
            class = "summary.pde_convergence")
}

#' @export
print.summary.pde_convergence <- function(x, ...) {
  cat("Mesh convergence study: h from ", fmt_num(x$h[2]), " to ", fmt_num(x$h[1]), ", ",
      x$ndof[1], " to ", x$ndof[2], " degrees of freedom\n", sep = "")
  print(x$rates, digits = 4, row.names = FALSE)
  invisible(x)
}

#' @export
print.pde_convergence <- function(x, ...) {
  cat("Mesh convergence study\n")
  print(as.data.frame(x), digits = 4, row.names = FALSE)
  invisible(x)
}

#' A posteriori error estimation
#'
#' `estimate_error()` computes the residual-based error estimator for a
#' linear elliptic solution:
#'
#' `eta_T^2 = h_T^2 ||f - c u_h - b . grad u_h + div(k grad u_h)||_T^2`
#' `+ 1/2 sum_{interior E in T} h_E ||[k grad u_h . n]||_E^2`
#' `+ sum_{Neumann/Robin E in T} h_E ||g - alpha u_h - k grad u_h . n||_E^2`
#' `+ sum_{Dirichlet E in T} h_E ||d/ds (g_D - u_h)||_E^2`.
#'
#' Edges with no condition are zero-flux Neumann edges. An edge under any
#' Dirichlet condition is a Dirichlet edge (as in the solver), and natural
#' conditions on the same edge add up. The Dirichlet term measures how well
#' the space represents the boundary data (tangential derivative by central
#' differences). The diffusion coefficient in the flux jumps is evaluated on
#' each side of an edge separately, so coefficients that jump across element
#' edges are handled; finite-element coefficients are allowed. For P1 with constant diffusion the divergence term vanishes;
#' variable scalar diffusion contributes `grad k . grad u_h`, computed with
#' central differences of `k` with a step of `1e-3` times the inradius of
#' each element, so the differences stay inside the element (pass
#' `diffusion_gradient` for an exact gradient). Region coefficients
#' ([region_coefficient()]) are constant on each element and contribute
#' nothing. The estimator is reliable and efficient up to
#' constants for diffusion-dominated problems (Verfurth 2013); it is not robust
#' for strongly advection-dominated problems.
#'
#' @param solution A `pde_solution` of an [elliptic_problem()].
#' @param diffusion_gradient Optional exact gradient of a scalar diffusion
#'   coefficient, `function(x, y)` returning a two-column matrix.
#' @return A list of class `pde_error_estimate` with `eta` (per element),
#'   `total`, and the interior and edge contributions.
#' @references Verfurth, R. (2013). *A Posteriori Error Estimation Techniques
#'   for Finite Element Methods*. Oxford University Press.
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 8)
#' sol <- solve_pde(poisson_problem(mesh, source = 1, boundary = dirichlet_bc(0)))
#' est <- estimate_error(sol)
#' est
#' @export
estimate_error <- function(solution, diffusion_gradient = NULL) {
  stopifnot(inherits(solution, "pde_solution"), inherits(solution$problem, "elliptic_problem"))
  prob <- solution$problem
  u <- solution$u; V <- u$space; mesh <- V$mesh
  geo <- element_geometry(V, 2L * V$degree + 2L)
  uh <- fe_qp(u, geo)
  R <- matrix(coef_qp(prob$source, geo, name = "source"), geo$m)
  cq <- coef_qp(prob$reaction, geo, name = "reaction")
  R <- R - cq * uh$value
  bq <- vec_qp(prob$advection, geo)
  if (!is.null(bq)) R <- R - (bq$x * uh$dx + bq$y * uh$dy)
  no_periodic(prob$bcs, "estimate_error()")
  k <- prob$diffusion
  if (inherits(k, "pde_tensor")) {
    stopf(paste0("estimate_error() supports scalar and constant tensor diffusion; for a spatially ",
                 "varying tensor use errornorm() against a reference solution."))
  }
  Kt <- if (is.matrix(k) && is.numeric(k)) k else NULL
  if (is.function(k) || inherits(k, "fe_function")) {
    kq <- coef_qp(k, geo, name = "diffusion")
    xx <- as.vector(geo$x); yy <- as.vector(geo$y)
    gk <- if (!is.null(diffusion_gradient)) {
      g <- as.matrix(diffusion_gradient(xx, yy))
      if (nrow(g) != length(xx) || ncol(g) != 2L || any(!is.finite(g))) {
        stopf("`diffusion_gradient` must return a finite two-column matrix with one row per point.")
      }
      g
    } else if (inherits(k, "pde_region_coefficient")) {
      matrix(0, length(xx), 2)  # constant per element
    } else if (is.function(k)) {
      # Element-local step: quadrature points lie well inside each element,
      # so a step of 1e-3 inradius never reaches a neighbor (or an interface
      # along element edges), whatever the length unit.
      P <- mesh$nodes; T <- mesh$triangles
      el <- function(a, b) sqrt(rowSums((P[T[, a], , drop = FALSE] - P[T[, b], , drop = FALSE])^2))
      inr <- abs(tri_det(P, T)) / (el(1, 2) + el(2, 3) + el(3, 1))
      hq <- rep(1e-3 * inr, geo$nq)
      kx <- function(x, y) eval_spatial(k, x, y, name = "diffusion")
      cbind((kx(xx + hq, yy) - kx(xx - hq, yy)) / (2 * hq),
            (kx(xx, yy + hq) - kx(xx, yy - hq)) / (2 * hq))
    } else {
      as.matrix(evaluate(k, xx, yy, gradient = TRUE)[, 2:3])
    }
    R <- R + matrix(gk[, 1], geo$m) * uh$dx + matrix(gk[, 2], geo$m) * uh$dy
  } else if (is.null(Kt)) {
    kq <- matrix(k, geo$m, geo$nq)
  }
  if (V$degree == 2L) {
    lap <- p2_div_flux(mesh, u, if (is.null(Kt)) diag(2) else Kt)
    R <- R + (if (is.null(Kt)) kq else 1) * lap
  }
  q <- mesh_quality(mesh)
  hT <- q$h_max
  eta_int <- hT^2 * as.vector((R^2 %*% geo$w) * geo$det)
  eta_edge <- edge_residuals(prob, u, Kt)
  eta <- sqrt(eta_int + eta_edge)
  structure(list(eta = eta, total = sqrt(sum(eta^2)), interior = sqrt(eta_int),
                 edges = sqrt(eta_edge)), class = "pde_error_estimate")
}

#' @export
print.pde_error_estimate <- function(x, ...) {
  cat("<pde_error_estimate> total ", fmt_num(x$total), " over ", length(x$eta),
      " elements (max ", fmt_num(max(x$eta)), ")\n", sep = "")
  invisible(x)
}

# div(K grad u_h) of a P2 function, constant per element (K constant).
p2_div_flux <- function(mesh, u, K) {
  P <- mesh$nodes; T <- mesh$triangles
  x1 <- P[T[, 1], 1]; y1 <- P[T[, 1], 2]
  b11 <- P[T[, 2], 1] - x1; b12 <- P[T[, 3], 1] - x1
  b21 <- P[T[, 2], 2] - y1; b22 <- P[T[, 3], 2] - y1
  det <- b11 * b22 - b12 * b21
  c11 <- b22 / det; c12 <- -b12 / det; c21 <- -b21 / det; c22 <- b11 / det
  # G = C K C^T
  g11 <- c11 * (K[1, 1] * c11 + K[1, 2] * c12) + c12 * (K[2, 1] * c11 + K[2, 2] * c12)
  g12 <- c11 * (K[1, 1] * c21 + K[1, 2] * c22) + c12 * (K[2, 1] * c21 + K[2, 2] * c22)
  g21 <- c21 * (K[1, 1] * c11 + K[1, 2] * c12) + c22 * (K[2, 1] * c11 + K[2, 2] * c12)
  g22 <- c21 * (K[1, 1] * c21 + K[1, 2] * c22) + c22 * (K[2, 1] * c21 + K[2, 2] * c22)
  Hxx <- c(4, 4, 0, -8, 0, 0); Hxy <- c(4, 0, 0, -4, 4, -4); Hyy <- c(4, 0, 4, 0, 0, -8)
  U <- matrix(u$values[u$space$cell_dofs], nrow(T))
  lap_phi <- outer(g11, Hxx) + outer(g12 + g21, Hxy) + outer(g22, Hyy)
  rowSums(lap_phi * U)
}

# Gradient of u_h at physical points (x, y) inside elements e.
grad_in_elements <- function(u, e, x, y) {
  mesh <- u$space$mesh; P <- mesh$nodes; T <- mesh$triangles[e, , drop = FALSE]
  x1 <- P[T[, 1], 1]; y1 <- P[T[, 1], 2]
  b11 <- P[T[, 2], 1] - x1; b12 <- P[T[, 3], 1] - x1
  b21 <- P[T[, 2], 2] - y1; b22 <- P[T[, 3], 2] - y1
  det <- b11 * b22 - b12 * b21
  dx <- x - x1; dy <- y - y1
  xi <- (b22 * dx - b12 * dy) / det; eta <- (-b21 * dx + b11 * dy) / det
  el <- u$space$element
  G <- el$grad(xi, eta)
  gr <- physical_gradients(mesh, e, G)
  U <- matrix(u$values[u$space$cell_dofs[e, , drop = FALSE]], length(e))
  list(dx = rowSums(gr$x * U), dy = rowSums(gr$y * U),
       value = rowSums(el$basis(xi, eta) * U))
}

edge_residuals <- function(prob, u, Kt) {
  mesh <- u$space$mesh; P <- mesh$nodes; E <- mesh$edges
  m <- nrow(mesh$triangles)
  gq <- quadrature_interval(u$space$degree + 2L)
  out <- numeric(m)
  ctr <- cbind((P[mesh$triangles[, 1], 1] + P[mesh$triangles[, 2], 1] + P[mesh$triangles[, 3], 1]) / 3,
               (P[mesh$triangles[, 1], 2] + P[mesh$triangles[, 2], 2] + P[mesh$triangles[, 3], 2]) / 3)
  # The diffusion coefficient is that of element `e`, so a coefficient that
  # jumps across an edge uses each side's own value: a region coefficient of
  # this mesh is read per element, any other coefficient is evaluated just
  # inside the element, at a point that stays distinct from the edge even at
  # large coordinates.
  kd <- prob$diffusion
  kel <- if (inherits(kd, "pde_region_coefficient") && identical(attr(kd, "mesh")$nodes, P) &&
             identical(attr(kd, "mesh")$triangles, mesh$triangles)) attr(kd, "element_values")
  flux_k <- function(x, y, g, e) {
    if (!is.null(Kt)) {
      list(x = Kt[1, 1] * g$dx + Kt[1, 2] * g$dy, y = Kt[2, 1] * g$dx + Kt[2, 2] * g$dy)
    } else {
      kv <- if (!is.null(kel)) kel[e] else {
        dx <- ctr[e, 1] - x; dy <- ctr[e, 2] - y
        fr <- pmin(0.5, pmax(1e-6, 16 * .Machine$double.eps * pmax(abs(x), abs(y)) / pmax(sqrt(dx^2 + dy^2), 1e-300)))
        eval_spatial(kd, x + fr * dx, y + fr * dy, name = "diffusion")
      }
      list(x = kv * g$dx, y = kv * g$dy)
    }
  }
  # interior edges
  int <- which(!is.na(mesh$edge_elements[, 2]))
  if (length(int)) {
    a <- E[int, 1]; b <- E[int, 2]
    len <- sqrt((P[b, 1] - P[a, 1])^2 + (P[b, 2] - P[a, 2])^2)
    nx <- (P[b, 2] - P[a, 2]) / len; ny <- -(P[b, 1] - P[a, 1]) / len
    t1 <- mesh$edge_elements[int, 1]; t2 <- mesh$edge_elements[int, 2]
    j2 <- numeric(length(int))
    for (k in seq_along(gq$weights)) {
      s <- gq$points[k, 1]
      x <- P[a, 1] + s * (P[b, 1] - P[a, 1]); y <- P[a, 2] + s * (P[b, 2] - P[a, 2])
      f1 <- flux_k(x, y, grad_in_elements(u, t1, x, y), t1)
      f2 <- flux_k(x, y, grad_in_elements(u, t2, x, y), t2)
      jump <- (f1$x - f2$x) * nx + (f1$y - f2$y) * ny
      j2 <- j2 + gq$weights[k] * len * jump^2
    }
    contrib <- len * j2
    out <- out + tabulate_sum(t1, 0.5 * contrib, m) + tabulate_sum(t2, 0.5 * contrib, m)
  }
  # boundary edges with natural conditions
  bd <- mesh$boundary
  if (nrow(bd)) {
    # As in the solver: an edge under any Dirichlet condition is a Dirichlet
    # edge; natural conditions on the same edge add up.
    is_dir <- rep(FALSE, nrow(bd))
    for (bc in prob$bcs) if (bc$type == "dirichlet") is_dir[boundary_rows(mesh, bc$marker)] <- TRUE
    sel <- which(!is_dir)
    if (length(sel)) {
      a <- bd$n1[sel]; b <- bd$n2[sel]; e <- bd$element[sel]
      nat <- Filter(function(bc) bc$type != "dirichlet", prob$bcs)
      rows <- lapply(nat, function(bc) sel %in% boundary_rows(mesh, bc$marker))
      len <- sqrt((P[b, 1] - P[a, 1])^2 + (P[b, 2] - P[a, 2])^2)
      nx <- (P[b, 2] - P[a, 2]) / len; ny <- -(P[b, 1] - P[a, 1]) / len
      r2 <- numeric(length(sel))
      for (k in seq_along(gq$weights)) {
        s <- gq$points[k, 1]
        x <- P[a, 1] + s * (P[b, 1] - P[a, 1]); y <- P[a, 2] + s * (P[b, 2] - P[a, 2])
        gi <- grad_in_elements(u, e, x, y)
        fl <- flux_k(x, y, gi, e)
        g <- numeric(length(sel)); al <- numeric(length(sel))
        for (i in seq_along(nat)) {
          bc <- nat[[i]]; w <- rows[[i]]
          if (!any(w)) next
          g[w] <- g[w] + eval_spatial(bc$value, x[w], y[w], name = "flux")
          if (bc$type == "robin") al[w] <- al[w] + eval_spatial(bc$alpha, x[w], y[w], name = "alpha")
        }
        res <- g - al * gi$value - (fl$x * nx + fl$y * ny)
        r2 <- r2 + gq$weights[k] * len * res^2
      }
      out <- out + tabulate_sum(e, len * r2, m)
    }
    # Dirichlet edges: approximation error of the boundary data,
    # h_E ||d/ds (g - u_h)||^2_E, so that data the space cannot represent are
    # not reported as zero error.
    dsel <- which(is_dir)
    if (length(dsel)) {
      a <- bd$n1[dsel]; b <- bd$n2[dsel]; e <- bd$element[dsel]
      tx <- P[b, 1] - P[a, 1]; ty <- P[b, 2] - P[a, 2]
      len <- sqrt(tx^2 + ty^2); tx <- tx / len; ty <- ty / len
      own <- rep(0L, length(dsel))
      for (i in seq_along(prob$bcs)) {
        if (prob$bcs[[i]]$type == "dirichlet") own[dsel %in% boundary_rows(mesh, prob$bcs[[i]]$marker)] <- i
      }
      d2 <- numeric(length(dsel))
      for (k in seq_along(gq$weights)) {
        s <- gq$points[k, 1]
        x <- P[a, 1] + s * tx * len; y <- P[a, 2] + s * ty * len
        gi <- grad_in_elements(u, e, x, y)
        dg <- numeric(length(dsel))
        # Central difference along the edge; the step is at least a few ulp
        # of the coordinates and the divisor is the distance actually taken.
        hs <- pmax(1e-6 * len, 8 * .Machine$double.eps * pmax(abs(x), abs(y)))
        xp <- x + hs * tx; yp <- y + hs * ty; xm <- x - hs * tx; ym <- y - hs * ty
        sep <- (xp - xm) * tx + (yp - ym) * ty
        for (i in unique(own)) {
          w <- own == i
          g <- prob$bcs[[i]]$value
          dg[w] <- (eval_spatial(g, xp[w], yp[w], name = "Dirichlet value") -
                      eval_spatial(g, xm[w], ym[w], name = "Dirichlet value")) / sep[w]
        }
        res <- dg - (gi$dx * tx + gi$dy * ty)
        d2 <- d2 + gq$weights[k] * len * res^2
      }
      out <- out + tabulate_sum(e, len * d2, m)
    }
  }
  out
}

tabulate_sum <- function(idx, val, n) {
  out <- numeric(n)
  s <- rowsum(val, idx, reorder = TRUE)
  out[as.integer(rownames(s))] <- s[, 1]
  out
}

#' Mark elements for refinement
#'
#' `mark_elements()` implements Doerfler (bulk) marking, which selects a
#' smallest set of elements whose squared indicators sum to at least
#' `theta` times the total, and maximum marking, which selects elements with
#' `eta_T >= theta * max(eta)`.
#'
#' @param eta Element error indicators.
#' @param theta Marking parameter in (0, 1].
#' @param strategy `"doerfler"` or `"maximum"`.
#' @return Integer indices of marked elements.
#' @references Doerfler, W. (1996). A convergent adaptive algorithm for
#'   Poisson's equation. *SIAM Journal on Numerical Analysis*, 33, 1106-1124.
#' @examples
#' mark_elements(c(0.1, 0.5, 0.2, 0.05), theta = 0.5)
#' @export
mark_elements <- function(eta, theta = 0.5, strategy = c("doerfler", "maximum")) {
  strategy <- match.arg(strategy)
  check_number(theta, "theta", positive = TRUE)
  if (inherits(eta, "pde_error_estimate")) eta <- eta$eta
  if (!is.numeric(eta) || !length(eta) || any(!is.finite(eta)) || any(eta < 0)) {
    stopf("`eta` must be a nonempty vector of finite nonnegative indicators.")
  }
  if (theta > 1) stopf("`theta` must be in (0, 1].")
  if (strategy == "maximum") return(which(eta >= theta * max(eta)))
  o <- order(eta, decreasing = TRUE)
  cs <- cumsum(eta[o]^2)
  k <- which(cs >= theta * sum(eta^2))[1]
  sort(o[seq_len(k)])
}

#' Adaptive finite-element solution
#'
#' Runs the loop solve, estimate, mark, refine for a linear elliptic problem,
#' using [estimate_error()], [mark_elements()] and newest vertex bisection
#' ([refine_mesh()]).
#'
#' @param problem An [elliptic_problem()].
#' @param theta Doerfler marking parameter.
#' @param max_iterations Maximum number of refinement steps.
#' @param max_dofs Stop when the number of degrees of freedom exceeds this.
#' @param tol Stop when the estimator falls below this.
#' @param exact,gradient Optional exact solution and gradient; when given,
#'   the true errors are recorded in the history.
#' @param ... Passed to [solve_pde()].
#' @return A list of class `pde_adaptive` with the final `solution`, the
#'   final `estimate` and a `history` data frame.
#' @examples
#' L <- pde_domain("polygon", vertices = cbind(c(0, 1, 1, 0.5, 0.5, 0), c(0, 0, 0.5, 0.5, 1, 1)))
#' prob <- poisson_problem(pde_mesh_polygon(L, h = 0.25), source = 1, boundary = dirichlet_bc(0))
#' res <- adaptive_solve(prob, max_iterations = 3)
#' res$history
#' @export
adaptive_solve <- function(problem, theta = 0.5, max_iterations = 8, max_dofs = 1e5,
                           tol = 0, exact = NULL, gradient = NULL, ...) {
  stopifnot(inherits(problem, "elliptic_problem"))
  hist <- list()
  prob <- problem
  for (it in 0:max_iterations) {
    sol <- solve_pde(prob, ...)
    est <- estimate_error(sol)
    row <- data.frame(iteration = it, ndof = sol$space$ndof,
                      elements = nrow(prob$mesh$triangles), estimator = est$total)
    if (!is.null(exact)) {
      row$error_L2 <- as.numeric(errornorm(sol, exact, "L2"))
      row$error_H1_semi <- as.numeric(errornorm(sol, exact, "H1_semi", gradient = gradient))
    }
    hist[[length(hist) + 1L]] <- row
    if (it == max_iterations || est$total <= tol || sol$space$ndof >= max_dofs) break
    marked <- mark_elements(est$eta, theta)
    prob <- update_mesh(prob, refine_mesh(prob$mesh, marked))
  }
  structure(list(solution = sol, estimate = est, history = do.call(rbind, hist)),
            class = "pde_adaptive")
}

#' @export
print.pde_adaptive <- function(x, ...) {
  cat("Adaptive finite-element solution\n")
  print(x$history, digits = 4, row.names = FALSE)
  invisible(x)
}
