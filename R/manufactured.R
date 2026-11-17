#' Manufactured solutions with symbolic derivatives
#'
#' `manufactured_solution()` takes an exact solution `u(x, y)` (or
#' `u(x, y, t)`) and coefficients written as strings or 'symcompR'
#' expressions, differentiates them symbolically with 'symcompR' and returns
#' vectorized R functions for the solution, its gradient and the source term
#'
#' `f = du/dt - div(k grad u) + b . grad u + c u`,
#'
#' where the time derivative is present only when `u` contains `t`. The
#' results plug directly into [elliptic_problem()], [heat_problem()],
#' [errornorm()] and [convergence_study()], which removes hand-derived source
#' terms as a source of error in verification studies.
#'
#' 'symcompR' is an optional dependency; the function stops with a message
#' when it is not installed.
#'
#' @param u Exact solution, a string such as `"sin(pi*x)*sin(pi*y)"` or a
#'   'symcompR' expression in `x`, `y` and optionally `t`.
#' @param diffusion Scalar diffusion coefficient `k(x, y)` (string, number or
#'   expression).
#' @param advection Optional advection field: a length-2 list or character
#'   vector of components.
#' @param reaction Reaction coefficient `c(x, y)`.
#' @return A list of class `pde_manufactured` with functions `exact`,
#'   `gradient` (returning a two-column matrix), `source`, `diffusion`,
#'   `advection` (or `NULL`) and `reaction`, the symbolic expressions in
#'   `symbolic`, and `time_dependent`. Time-dependent functions take
#'   `(x, y, t)`.
#' @examples
#' if (requireNamespace("symcompR", quietly = TRUE)) {
#'   ms <- manufactured_solution("sin(pi*x)*sin(pi*y)", diffusion = "1 + x")
#'   prob <- elliptic_problem(pde_mesh_rectangle(nx = 8), diffusion = ms$diffusion,
#'                            source = ms$source, bcs = dirichlet_bc(0))
#'   errornorm(solve_pde(prob), ms$exact, "H1_semi", gradient = ms$gradient)
#' }
#' @export
manufactured_solution <- function(u, diffusion = 1, advection = NULL, reaction = 0) {
  if (!requireNamespace("symcompR", quietly = TRUE)) {
    stopf("`manufactured_solution()` needs the 'symcompR' package; install it first.")
  }
  sx <- symcompR::Symbol("x"); sy <- symcompR::Symbol("y"); st <- symcompR::Symbol("t")
  as_sym <- function(e) symcompR::sympify(if (is.numeric(e)) as.character(e) else e)
  U <- as_sym(u); K <- as_sym(diffusion); C <- as_sym(reaction)
  B <- NULL
  if (!is.null(advection)) {
    if (length(advection) != 2L) stopf("`advection` must have two components.")
    B <- lapply(advection, as_sym)
  }
  syms <- function(e) vapply(symcompR::free_symbols(e), format, character(1))
  for (e in c(list(U, K, C), B)) {
    bad <- setdiff(syms(e), c("x", "y", "t"))
    if (length(bad)) stopf("Unknown symbol(s) %s; expressions may use only x, y and t.", paste(bad, collapse = ", "))
  }
  td <- "t" %in% syms(U)
  for (e in c(list(K, C), B)) {
    if ("t" %in% syms(e)) stopf("Coefficients (diffusion, advection, reaction) must not depend on time.")
  }
  ux <- diff(U, sx); uy <- diff(U, sy)
  f <- -(diff(K * ux, sx) + diff(K * uy, sy)) + C * U
  if (!is.null(B)) f <- f + B[[1]] * ux + B[[2]] * uy
  if (td) f <- f + diff(U, st)
  f <- symcompR::simplify(f)
  args <- if (td) list(sx, sy, st) else list(sx, sy)
  fun <- function(e, time = td) {
    g <- symcompR::lambdify(if (time) args else list(sx, sy), e)
    if (time) {
      function(x, y, t) rep_len(as.numeric(g(x, y, t)), length(x))
    } else {
      function(x, y) rep_len(as.numeric(g(x, y)), length(x))
    }
  }
  gx <- fun(ux); gy <- fun(uy)
  gradient <- if (td) function(x, y, t) cbind(gx(x, y, t), gy(x, y, t))
              else function(x, y) cbind(gx(x, y), gy(x, y))
  adv <- NULL
  if (!is.null(B)) {
    bx <- fun(B[[1]], FALSE); by <- fun(B[[2]], FALSE)
    adv <- function(x, y) cbind(bx(x, y), by(x, y))
  }
  structure(list(exact = fun(U), gradient = gradient, source = fun(f),
                 diffusion = fun(K, FALSE), advection = adv, reaction = fun(C, FALSE),
                 symbolic = list(u = U, source = f, diffusion = K, advection = B, reaction = C),
                 time_dependent = td),
            class = "pde_manufactured")
}

#' @export
print.pde_manufactured <- function(x, ...) {
  cat("<pde_manufactured>", if (x$time_dependent) "time-dependent" else "stationary", "\n")
  cat("  u = ", format(x$symbolic$u), "\n", sep = "")
  cat("  f = ", format(x$symbolic$source), "\n", sep = "")
  invisible(x)
}
