#' Eigenvalue problems
#'
#' Computes the eigenpairs `K x = lambda M x` of the discrete operator of an
#' [elliptic_problem()] (diffusion and reaction, with Dirichlet conditions
#' removing the prescribed degrees of freedom; Neumann and Robin terms enter
#' `K`) or an [elasticity_problem()] (vibration modes, with `M` the mass
#' matrix of unit density), nearest to the shift `sigma`. Periodic
#' conditions reduce both matrices as in [solve_pde()], and the returned
#' eigenvectors satisfy them.
#'
#' The method is shift-and-invert subspace iteration with Rayleigh-Ritz
#' projection: one sparse factorization of `K - sigma M`, a block of
#' `n + 8` vectors solved with that factorization at each iteration, and a
#' dense generalized eigenproblem of the block size. It stops when every
#' requested pair has normwise relative residual
#' `||K x - lambda M x|| / ((||K|| + |lambda| ||M||) ||x||)` below `tol`. Problems without Dirichlet
#' conditions have a zero eigenvalue; the default shift is then slightly
#' negative so that `K - sigma M` is nonsingular.
#'
#' @param problem An [elliptic_problem()] without advection, or an
#'   [elasticity_problem()]. Sources, Neumann data and body forces are
#'   ignored; only the operator matters.
#' @param n Number of eigenpairs.
#' @param sigma Shift; the eigenvalues nearest to it are returned.
#' @param tol Relative residual tolerance.
#' @param maxit Maximum number of iterations.
#' @return A list of class `pde_eigen` with `values` (ascending distance to
#'   `sigma`), `vectors` (a list of [fe_function()]s, or of displacement
#'   pairs for elasticity, normalized to `x'Mx = 1`), `residuals` and
#'   `iterations`.
#' @examples
#' p <- elliptic_problem(pde_mesh_rectangle(nx = 16), bcs = dirichlet_bc(0))
#' ev <- solve_eigen(p, n = 4)
#' ev$values / pi^2   # close to 2, 5, 5, 8
#' @export
solve_eigen <- function(problem, n = 6, sigma = NULL, tol = 1e-10, maxit = 500) {
  n <- check_count(n, "n", min = 1L)
  maxit <- check_count(maxit, "maxit", min = 1L)
  check_tol(tol, "tol")
  if (inherits(problem, "elliptic_problem")) {
    if (!is.null(problem$advection)) stopf("Eigenproblems need a symmetric operator: remove the advection.")
    if (problem$stabilization != "none") stopf("Eigenproblems do not support stabilization.")
    disc <- discretize(problem)
    M <- assemble_matrix(disc$space, "mass")
  } else if (inherits(problem, "elasticity_problem")) {
    if (identical(problem$formulation, "mixed")) {
      stopf("Eigenproblems use the displacement formulation; set `formulation = \"displacement\"`.")
    }
    disc <- discretize(problem)
    M1 <- assemble_matrix(disc$space, "mass")
    nd <- disc$space$ndof
    Z <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(nd, nd))
    M <- rbind(cbind(M1, Z), cbind(Z, M1))
  } else {
    stopf("`problem` must be an elliptic_problem() or an elasticity_problem().")
  }
  A <- disc$A; f <- disc$free; Pm <- disc$periodic
  if (!is.null(Pm)) {
    # Periodic constraints: the same reduction as solve_pde(), applied to
    # both matrices; eigenvectors are lifted back with P.
    cs <- constrained_system(disc)
    A <- cs$A; f <- cs$free
    M <- Matrix::drop0(Matrix::crossprod(Pm, M %*% Pm))
  }
  K <- A[f, f, drop = FALSE]; Mf <- M[f, f, drop = FALSE]
  if (!scaled_symmetric(K, tol = 1e-10)) stopf("Eigenproblems need a symmetric operator.")
  N <- length(f)
  # With n close to N the block spans the whole space and Rayleigh-Ritz is
  # exact, so the full spectrum (n == N) is allowed.
  if (n > N) stopf("At most %d eigenpairs exist for %d free degrees of freedom.", N, N)
  if (is.null(sigma)) {
    sigma <- if (length(disc$dirichlet$dofs)) 0 else -1e-3 * mean(Matrix::diag(K)) / mean(Matrix::diag(Mf))
  }
  check_number(sigma, "sigma")
  slv <- linear_solver(K - sigma * Mf)
  normK <- max(Matrix::rowSums(abs(K))); normM <- max(Matrix::rowSums(abs(Mf)))
  p <- min(N, n + 8L)
  # Reproducible random start block; the user's random number stream is
  # restored afterwards.
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(20261117L)
  X <- matrix(stats::runif(N * p) - 0.5, N, p)
  if (had) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  it <- 0L
  repeat {
    it <- it + 1L
    S <- slv$solve(as.matrix(Mf %*% X))
    if (!solve_ok(S$info)) warn_solver(S$info, "shifted solve")
    Y <- qr.Q(qr(S$x))
    Kr <- as.matrix(Matrix::crossprod(Y, K %*% Y)); Mr <- as.matrix(Matrix::crossprod(Y, Mf %*% Y))
    Kr <- (Kr + t(Kr)) / 2; Mr <- (Mr + t(Mr)) / 2
    R <- chol(Mr)
    Ri <- backsolve(R, diag(p))
    e <- eigen(t(Ri) %*% Kr %*% Ri, symmetric = TRUE)
    o <- order(abs(e$values - sigma))
    lam <- e$values[o]
    X <- Y %*% (Ri %*% e$vectors[, o, drop = FALSE])
    KX <- as.matrix(K %*% X[, seq_len(n), drop = FALSE])
    xn <- sqrt(colSums(X[, seq_len(n), drop = FALSE]^2))
    res <- sqrt(colSums((KX - as.matrix(Mf %*% X[, seq_len(n), drop = FALSE]) %*% diag(lam[seq_len(n)], n))^2)) /
      pmax((normK + abs(lam[seq_len(n)]) * normM) * xn, 1e-300)
    if (all(res <= tol) || it >= maxit) break
  }
  if (any(res > tol)) warnf("Eigenpairs did not converge in %d iterations (largest residual %s).", it, fmt_num(max(res)))
  full <- function(x) {
    v <- numeric(nrow(M)); v[f] <- x
    if (is.null(Pm)) v else as.numeric(Pm %*% v)
  }
  vecs <- lapply(seq_len(n), function(j) {
    v <- full(X[, j])
    if (inherits(problem, "elasticity_problem")) {
      nd <- disc$space$ndof
      list(ux = fe_function(disc$space, v[seq_len(nd)], "ux"), uy = fe_function(disc$space, v[nd + seq_len(nd)], "uy"))
    } else fe_function(disc$space, v, sprintf("mode %d", j))
  })
  structure(list(values = lam[seq_len(n)], vectors = vecs, residuals = res, iterations = it, sigma = sigma),
            class = "pde_eigen")
}

#' @export
print.pde_eigen <- function(x, ...) {
  cat("<pde_eigen> ", length(x$values), " eigenpairs nearest ", fmt_num(x$sigma), " (", x$iterations,
      " iterations, largest residual ", fmt_num(max(x$residuals), 2), ")\n", sep = "")
  cat("  values:", fmt_num(x$values, 8), "\n")
  invisible(x)
}
