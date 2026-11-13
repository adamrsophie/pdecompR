#' Weak-form language
#'
#' A small language for linear scalar variational problems. `trial_function()`
#' and `test_function()` create the symbols `u` and `v`; `grad()` and `dot()`
#' build gradients and inner products; integrands are combined with `*`,
#' `+` and `-`, and turned into forms by `integrate_domain()` (`dx`) and
#' `integrate_boundary()` (`ds` over marked boundary edges).
#'
#' Supported integrands, with `c` a scalar coefficient (number or
#' `function(x, y)`) and `b` a vector coefficient (length-2 vector or
#' `function(x, y)` returning a two-column matrix):
#'
#' * bilinear, domain: `c * dot(grad(u), grad(v))`, `c * u * v`,
#'   `dot(b, grad(u)) * v`, `u * dot(b, grad(v))`;
#' * bilinear, boundary: `c * u * v`;
#' * linear, domain: `c * v`, `dot(b, grad(v))`;
#' * linear, boundary: `c * v`.
#'
#' Scalar coefficients multiply; products are evaluated pointwise at the
#' quadrature points. Write a coefficient function first and multiply it
#' afterwards (`f <- function(x, y) x * y; f * v`): in
#' `function(x, y) x * y * v` the `* v` belongs to the function body. The
#' functions `grad()` and `dot()` only act on these
#' symbolic objects and do not mask any base R function.
#'
#' @param space An [fe_space()].
#' @param x,a,b Weak-form expressions or coefficients.
#' @param expr A weak-form integrand.
#' @param marker Boundary markers for `integrate_boundary()`.
#' @return Expressions of class `fe_expr`, forms of class `fe_form`.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 8))
#' u <- trial_function(V); v <- test_function(V)
#' a <- integrate_domain(dot(grad(u), grad(v)))
#' f <- function(x, y) 2 * pi^2 * sin(pi * x) * sin(pi * y)
#' L <- integrate_domain(f * v)
#' prob <- variational_problem(a, L, dirichlet_bc(0))
#' sol <- solve_pde(prob)
#' @export
trial_function <- function(space) {
  stopifnot(inherits(space, "fe_space"))
  new_expr(list(mono(trial = "value")), space)
}

#' @rdname trial_function
#' @export
test_function <- function(space) {
  stopifnot(inherits(space, "fe_space"))
  new_expr(list(mono(test = "value")), space)
}

mono <- function(trial = "none", test = "none", coefs = list(), vec = NULL, rank = 0L,
                 contract = "none") {
  list(trial = trial, test = test, coefs = coefs, vec = vec, rank = rank, contract = contract)
}

new_expr <- function(terms, space) structure(list(terms = terms, space = space), class = "fe_expr")

as_expr <- function(x, space) {
  if (inherits(x, "fe_expr")) return(x)
  if ((is.numeric(x) && length(x) == 1L) || is.function(x) || inherits(x, "fe_function")) {
    return(new_expr(list(mono(coefs = list(x))), space))
  }
  if (is.numeric(x) && length(x) == 2L) {
    return(new_expr(list(mono(vec = x, rank = 1L)), space))
  }
  stopf("Cannot use an object of class '%s' in a weak form.", class(x)[1])
}

#' @rdname trial_function
#' @export
grad <- function(x) {
  if (!inherits(x, "fe_expr")) stopf("`grad()` applies to trial or test functions.")
  terms <- lapply(x$terms, function(t) {
    if (t$rank != 0L || !is.null(t$vec)) stopf("`grad()` applies to u or v only.")
    if (t$trial == "value" && t$test == "none") t$trial <- "grad"
    else if (t$test == "value" && t$trial == "none") t$test <- "grad"
    else stopf("`grad()` applies to u or v only.")
    t$rank <- 1L
    t
  })
  new_expr(terms, x$space)
}

#' @rdname trial_function
#' @export
vector_coefficient <- function(x) {
  if (!(is.function(x) || (is.numeric(x) && length(x) == 2L))) {
    stopf("A vector coefficient is a length-2 vector or a function returning two columns.")
  }
  structure(list(value = x), class = "fe_vector_coefficient")
}

#' @rdname trial_function
#' @export
dot <- function(a, b) {
  sp <- if (inherits(a, "fe_expr")) a$space else if (inherits(b, "fe_expr")) b$space else NULL
  if (is.null(sp)) stopf("`dot()` needs a weak-form expression.")
  vecify <- function(z) {
    if (inherits(z, "fe_vector_coefficient")) return(new_expr(list(mono(vec = z$value, rank = 1L)), sp))
    as_expr(z, sp)
  }
  A <- vecify(a); B <- vecify(b)
  terms <- list()
  for (s in A$terms) for (t in B$terms) {
    if (s$rank != 1L || t$rank != 1L) stopf("`dot()` needs two vector-valued arguments.")
    r <- mono(trial = if (s$trial != "none") s$trial else t$trial,
              test = if (s$test != "none") s$test else t$test,
              coefs = c(s$coefs, t$coefs))
    if (s$trial != "none" && t$trial != "none") stopf("A term contains the trial function twice.")
    if (s$test != "none" && t$test != "none") stopf("A term contains the test function twice.")
    vec <- if (!is.null(s$vec)) s$vec else t$vec
    if (!is.null(s$vec) && !is.null(t$vec)) stopf("`dot()` of two coefficients is not supported.")
    if (r$trial == "grad" && r$test == "grad") r$contract <- "gg"
    else if (r$trial == "grad") { r$contract <- "bu"; r$vec <- vec }
    else if (r$test == "grad") { r$contract <- "bv"; r$vec <- vec }
    else stopf("Unsupported `dot()` arguments.")
    terms[[length(terms) + 1L]] <- r
  }
  new_expr(terms, sp)
}

#' @export
Ops.fe_expr <- function(e1, e2) {
  sp <- if (inherits(e1, "fe_expr")) e1$space else e2$space
  if (missing(e2)) {
    if (.Generic == "-") return(as_expr(-1, sp) * e1)
    if (.Generic == "+") return(e1)
  }
  switch(.Generic,
    "+" = new_expr(c(as_expr(e1, sp)$terms, as_expr(e2, sp)$terms), sp),
    "-" = new_expr(c(as_expr(e1, sp)$terms, as_expr(-1 * as_expr(e2, sp), sp)$terms), sp),
    "*" = {
      A <- as_expr(e1, sp); B <- as_expr(e2, sp)
      terms <- list()
      for (s in A$terms) for (t in B$terms) {
        if (s$trial != "none" && t$trial != "none") stopf("A term contains the trial function twice.")
        if (s$test != "none" && t$test != "none") stopf("A term contains the test function twice.")
        if (s$rank + t$rank > 1L) stopf("Use `dot()` to multiply vector-valued expressions.")
        r <- mono(trial = if (s$trial != "none") s$trial else t$trial,
                  test = if (s$test != "none") s$test else t$test,
                  coefs = c(s$coefs, t$coefs),
                  vec = if (!is.null(s$vec)) s$vec else t$vec,
                  rank = s$rank + t$rank,
                  contract = if (s$contract != "none") s$contract else t$contract)
        terms[[length(terms) + 1L]] <- r
      }
      new_expr(terms, sp)
    },
    stopf("Operator `%s` is not defined for weak-form expressions.", .Generic))
}

#' @export
print.fe_expr <- function(x, ...) {
  cat("<fe_expr> with", length(x$terms), "term(s)\n")
  invisible(x)
}

#' @rdname trial_function
#' @export
integrate_domain <- function(expr) {
  if (!inherits(expr, "fe_expr")) stopf("`expr` must be a weak-form expression.")
  structure(list(parts = list(list(where = "domain", marker = NULL, expr = expr)),
                 space = expr$space), class = "fe_form")
}

#' @rdname trial_function
#' @export
integrate_boundary <- function(expr, marker = "boundary") {
  if (!inherits(expr, "fe_expr")) stopf("`expr` must be a weak-form expression.")
  boundary_rows(expr$space$mesh, marker)
  structure(list(parts = list(list(where = "boundary", marker = marker, expr = expr)),
                 space = expr$space), class = "fe_form")
}

#' @export
Ops.fe_form <- function(e1, e2) {
  if (.Generic == "+" && inherits(e1, "fe_form") && inherits(e2, "fe_form")) {
    return(structure(list(parts = c(e1$parts, e2$parts), space = e1$space), class = "fe_form"))
  }
  if (.Generic == "-" && inherits(e1, "fe_form") && inherits(e2, "fe_form")) {
    neg <- lapply(e2$parts, function(p) { p$expr <- -1 * p$expr; p })
    return(structure(list(parts = c(e1$parts, neg), space = e1$space), class = "fe_form"))
  }
  stopf("Forms can only be added or subtracted.")
}

#' @export
print.fe_form <- function(x, ...) {
  cat("<fe_form> with", length(x$parts), "integral(s)\n")
  invisible(x)
}

# Product of scalar coefficients as one function of (x, y) or NULL.
coef_product <- function(coefs) {
  if (!length(coefs)) return(1)
  if (all(vapply(coefs, function(z) is.numeric(z), logical(1)))) return(prod(unlist(coefs)))
  function(x, y) {
    out <- rep(1, length(x))
    for (z in coefs) {
      out <- out * if (inherits(z, "fe_function")) evaluate(z, x, y) else eval_spatial(z, x, y)
    }
    out
  }
}

vec_scale <- function(vec, s) {
  if (identical(s, 1)) return(vec)
  function(x, y) {
    v <- if (is.numeric(vec)) matrix(vec, length(x), 2, byrow = TRUE) else as.matrix(vec(x, y))
    v * eval_spatial(s, x, y)
  }
}

# Assemble the matrix or vector of a form.
assemble_form <- function(form, kind = c("bilinear", "linear")) {
  kind <- match.arg(kind)
  V <- form$space
  n <- V$ndof
  out <- if (kind == "bilinear") Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                                                       dims = c(n, n)) else numeric(n)
  geo <- NULL
  for (p in form$parts) {
    for (t in p$expr$terms) {
      bil <- t$trial != "none"
      if (t$test == "none") stopf("Every term of a form must contain the test function.")
      if (bil != (kind == "bilinear")) {
        stopf("A %s form contains a %s term.", kind, if (bil) "bilinear" else "linear")
      }
      if (t$rank != 0L) stopf("A form integrand must be scalar; use `dot()`.")
      s <- coef_product(t$coefs)
      if (p$where == "boundary") {
        if (t$trial == "grad" || t$test == "grad") stopf("Boundary integrands may not contain gradients.")
        if (bil) out <- out + assemble_boundary_matrix(V, s, p$marker)
        else out <- out + assemble_boundary_vector(V, s, p$marker)
        next
      }
      if (is.null(geo)) geo <- element_geometry(V)
      if (bil) {
        A <- switch(paste(t$trial, t$test),
          "grad grad" = local_bilinear(geo, diff = coef_qp(s, geo)),
          "value value" = local_bilinear(geo, mass = coef_qp(s, geo)),
          "grad value" = local_bilinear(geo, adv = vec_qp(vec_scale(t$vec, s), geo)),
          "value grad" = local_bilinear(geo, advt = vec_qp(vec_scale(t$vec, s), geo)),
          stopf("Unsupported bilinear term."))
        out <- out + global_matrix(V, A)
      } else {
        b <- if (t$test == "value") local_linear(geo, f = coef_qp(s, geo))
             else local_linear(geo, g = vec_qp(vec_scale(t$vec, s), geo))
        out <- out + global_vector(V, b)
      }
    }
  }
  out
}

#' Variational problems
#'
#' `variational_problem()` pairs a bilinear form `a(u, v)` and a linear form
#' `L(v)` from the weak-form language with Dirichlet conditions (Neumann and
#' Robin terms are written directly as boundary integrals). Solving it
#' assembles both forms, eliminates the Dirichlet degrees of freedom and
#' solves the linear system.
#'
#' @param bilinear An `fe_form` containing the trial and test functions.
#' @param linear An `fe_form` containing only the test function.
#' @param bcs Dirichlet conditions ([dirichlet_bc()]).
#' @return An object of class `variational_problem`.
#' @export
variational_problem <- function(bilinear, linear, bcs = list()) {
  stopifnot(inherits(bilinear, "fe_form"), inherits(linear, "fe_form"))
  bcs <- normalize_bcs(bcs)
  if (any(vapply(bcs, function(b) b$type != "dirichlet", logical(1)))) {
    stopf("Write Neumann and Robin terms as boundary integrals; only Dirichlet conditions are allowed here.")
  }
  check_bc_markers(bcs, bilinear$space$mesh)
  structure(list(a = bilinear, L = linear, bcs = bcs, space = bilinear$space,
                 mesh = bilinear$space$mesh, degree = bilinear$space$degree),
            class = c("variational_problem", "pde_problem"))
}

#' @export
discretize.variational_problem <- function(problem, ...) {
  t0 <- proc.time()[["elapsed"]]
  V <- problem$space
  A <- assemble_form(problem$a, "bilinear")
  b <- assemble_form(problem$L, "linear")
  dd <- dirichlet_data(V, problem$bcs)
  structure(list(space = V, A = A, b = b, dirichlet = dd,
                 free = setdiff(seq_len(V$ndof), dd$dofs), pure_neumann = FALSE,
                 assembly_time = proc.time()[["elapsed"]] - t0),
            class = "pde_discretization")
}

#' @export
solve_pde.variational_problem <- function(problem, method = "auto", ...) {
  sol <- solve_discretization(discretize(problem), method = method, ...)
  sol$problem <- problem
  sol
}
