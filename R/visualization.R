#' Plot meshes and finite-element functions
#'
#' Base-graphics plots that work in any R session (the solvers themselves
#' never need a graphics device).
#'
#' * `plot.pde_mesh()` draws the triangles and colors boundary edges by
#'   marker.
#' * `plot.fe_function()` draws a filled plot (each triangle colored by the
#'   mean of its vertex values) or contour lines from a regular-grid
#'   evaluation (`type = "contour"`).
#' * `plot.pde_solution()` plots the solution; `plot.pde_transient()` plots a
#'   stored time; `plot.pde_convergence()` draws errors against `h` on log
#'   scales with reference slopes; `plot_indicators()` draws element-wise
#'   values such as error indicators.
#'
#' @param x Object to plot.
#' @param markers Color boundary edges by marker.
#' @param type `"fill"` or `"contour"`.
#' @param col Color palette (vector of colors).
#' @param nlevels Number of contour levels or color levels.
#' @param mesh Draw mesh edges on top of a filled plot.
#' @param time Time to plot for transient solutions.
#' @param species Species for reaction-diffusion solutions.
#' @param main Title.
#' @param values Element-wise values for `plot_indicators()`.
#' @param ... Passed to the underlying plotting functions.
#' @return The input, invisibly.
#' @name pde_plots
#' @examples
#' mesh <- pde_mesh_rectangle(nx = 6)
#' plot(mesh)
#' u <- interpolate(function(x, y) sin(pi * x) * y, fe_space(mesh))
#' plot(u)
#' plot(u, type = "contour")
NULL

#' @rdname pde_plots
#' @export
plot.pde_mesh <- function(x, markers = TRUE, main = "Mesh", ...) {
  P <- x$nodes; T <- x$triangles
  plot(range(P[, 1]), range(P[, 2]), type = "n", asp = 1, xlab = "x", ylab = "y",
       main = main, ...)
  E <- x$edges
  segments(P[E[, 1], 1], P[E[, 1], 2], P[E[, 2], 1], P[E[, 2], 2], col = "grey55", lwd = 0.6)
  if (markers) {
    bd <- x$boundary
    mk <- unique(bd$marker)
    cols <- grDevices::hcl.colors(max(2, length(mk)), "Dark 3")[seq_along(mk)]
    segments(P[bd$n1, 1], P[bd$n1, 2], P[bd$n2, 1], P[bd$n2, 2],
             col = cols[match(bd$marker, mk)], lwd = 2.5)
    legend("topright", legend = mk, col = cols, lwd = 2.5, bty = "n", cex = 0.8, inset = c(0, -0.02))
  }
  invisible(x)
}

#' @rdname pde_plots
#' @export
plot.fe_function <- function(x, type = c("fill", "contour"), col = grDevices::hcl.colors(64, "viridis"),
                             nlevels = 12, mesh = FALSE, main = x$name, ...) {
  type <- match.arg(type)
  msh <- x$space$mesh
  P <- msh$nodes; T <- msh$triangles
  nv <- nrow(P)
  vv <- x$values[seq_len(nv)]
  if (type == "fill") {
    tv <- (vv[T[, 1]] + vv[T[, 2]] + vv[T[, 3]]) / 3
    rng <- range(tv)
    k <- if (diff(rng) > 0) pmin(length(col), 1L + floor((tv - rng[1]) / diff(rng) * length(col))) else rep(1L, length(tv))
    op <- par(no.readonly = TRUE); on.exit(par(op))
    graphics::layout(matrix(1:2, 1), widths = c(5, 1))
    par(mar = c(4, 4, 3, 1))
    plot(range(P[, 1]), range(P[, 2]), type = "n", asp = 1, xlab = "x", ylab = "y", main = main, ...)
    px <- rbind(P[T[, 1], 1], P[T[, 2], 1], P[T[, 3], 1], NA)
    py <- rbind(P[T[, 1], 2], P[T[, 2], 2], P[T[, 3], 2], NA)
    polygon(as.vector(px), as.vector(py), col = col[k], border = if (mesh) "grey40" else col[k],
            lwd = if (mesh) 0.4 else 0.2)
    par(mar = c(4, 0.5, 3, 3.5))
    z <- seq(rng[1], rng[2], length.out = length(col))
    image(1, z, matrix(z, 1), col = col, axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1, cex.axis = 0.8)
    box()
  } else {
    n <- 101L
    gx <- seq(min(P[, 1]), max(P[, 1]), length.out = n)
    gy <- seq(min(P[, 2]), max(P[, 2]), length.out = n)
    g <- expand.grid(x = gx, y = gy)
    z <- matrix(evaluate(x, g$x, g$y), n, n)
    plot(range(P[, 1]), range(P[, 2]), type = "n", asp = 1, xlab = "x", ylab = "y", main = main, ...)
    bd <- msh$boundary
    segments(P[bd$n1, 1], P[bd$n1, 2], P[bd$n2, 1], P[bd$n2, 2], col = "grey30")
    contour(gx, gy, z, nlevels = nlevels, add = TRUE, col = col[round(seq(1, length(col), length.out = nlevels))])
  }
  invisible(x)
}

#' @rdname pde_plots
#' @export
plot.pde_solution <- function(x, ...) {
  plot(x$u, ...)
  invisible(x)
}

#' @rdname pde_plots
#' @export
plot.pde_transient <- function(x, time = max(x$times), species = 1, ...) {
  u <- solution_at(x, time, species)
  plot(u, main = sprintf("%s at t = %s", u$name, fmt_num(time)), ...)
  invisible(x)
}

#' @rdname pde_plots
#' @export
plot.pde_convergence <- function(x, main = "Mesh convergence", ...) {
  norms <- attr(x, "norms")
  ys <- unlist(x[norms])
  plot(x$h, x[[norms[1]]], log = "xy", type = "n", ylim = range(ys), xlab = "h",
       ylab = "error", main = main, ...)
  cols <- grDevices::hcl.colors(max(2, length(norms)), "Dark 3")
  for (i in seq_along(norms)) {
    lines(x$h, x[[norms[i]]], type = "b", pch = 19, col = cols[i])
    r <- round(utils::tail(x[[paste0("rate_", norms[i])]], 1))
    if (is.finite(r) && r > 0) {
      ref <- x[[norms[i]]][1] * (x$h / x$h[1])^r
      lines(x$h, ref, lty = 3, col = cols[i])
    }
  }
  legend("bottomright", legend = norms, col = cols[seq_along(norms)], pch = 19, lty = 1, bty = "n")
  invisible(x)
}

#' @rdname pde_plots
#' @export
plot_indicators <- function(mesh, values, col = grDevices::hcl.colors(64, "inferno", rev = TRUE),
                            main = "Element indicators", ...) {
  if (inherits(values, "pde_error_estimate")) values <- values$eta
  stopifnot(inherits(mesh, "pde_mesh"), length(values) == nrow(mesh$triangles))
  P <- mesh$nodes; T <- mesh$triangles
  rng <- range(values)
  k <- if (diff(rng) > 0) pmin(length(col), 1L + floor((values - rng[1]) / diff(rng) * length(col))) else rep(1L, length(values))
  plot(range(P[, 1]), range(P[, 2]), type = "n", asp = 1, xlab = "x", ylab = "y", main = main, ...)
  px <- rbind(P[T[, 1], 1], P[T[, 2], 1], P[T[, 3], 1], NA)
  py <- rbind(P[T[, 1], 2], P[T[, 2], 2], P[T[, 3], 2], NA)
  polygon(as.vector(px), as.vector(py), col = col[k], border = "grey60", lwd = 0.3)
  invisible(mesh)
}

#' @rdname pde_plots
#' @export
plot.pde_adaptive <- function(x, main = "Adaptive refinement", ...) {
  h <- x$history
  cols <- c("black", grDevices::hcl.colors(3, "Dark 3"))
  ys <- unlist(h[intersect(c("estimator", "error_H1_semi", "error_L2"), names(h))])
  plot(h$ndof, h$estimator, log = "xy", type = "b", pch = 19, ylim = range(ys),
       xlab = "degrees of freedom", ylab = "error / estimator", main = main, ...)
  leg <- "estimator"
  if (!is.null(h$error_H1_semi)) {
    lines(h$ndof, h$error_H1_semi, type = "b", pch = 17, col = cols[2]); leg <- c(leg, "H1 error")
  }
  if (!is.null(h$error_L2)) {
    lines(h$ndof, h$error_L2, type = "b", pch = 15, col = cols[3]); leg <- c(leg, "L2 error")
  }
  legend("topright", legend = leg, col = cols[seq_along(leg)], pch = c(19, 17, 15)[seq_along(leg)],
         lty = 1, bty = "n")
  invisible(x)
}
