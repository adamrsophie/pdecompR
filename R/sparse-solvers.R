#' Sparse linear solvers
#'
#' `solve_linear_system()` solves `A x = b` for a sparse matrix and reports
#' diagnostics. `linear_solver()` prepares a reusable solver (for example a
#' Cholesky factorization) for repeated right-hand sides, as in time stepping.
#'
#' Methods:
#'
#' * `"auto"`: sparse Cholesky (CHOLMOD through `Matrix`) for symmetric
#'   matrices, falling back to sparse LU when the factorization shows the
#'   matrix is not positive definite; sparse LU for nonsymmetric matrices.
#' * `"cholesky"`: sparse Cholesky; fails for matrices that are not
#'   symmetric positive definite.
#' * `"lu"`: sparse LU (any nonsingular matrix).
#' * `"cg"`: preconditioned conjugate gradients (symmetric positive definite
#'   matrices only).
#' * `"bicgstab"`, `"gmres"`: preconditioned BiCGSTAB and restarted GMRES for
#'   general matrices.
#'
#' Iterative methods use Jacobi (diagonal) preconditioning by default and stop
#' when `||b - A x|| <= tol * ||b||`. Direct methods report the same relative
#' residual, and status `"singular"` when the ratio of the smallest to the
#' largest pivot of the factorization is below about `2e-13` (the matrix is
#' singular to working precision, for example a stiffness matrix without
#' boundary conditions). A small residual says the linear system was solved accurately;
#' it says nothing about the discretization error of the PDE.
#'
#' @param A A square sparse matrix (`Matrix` class) or base matrix.
#' @param b Right-hand side vector.
#' @param method Solver, see Details.
#' @param tol Relative residual tolerance for iterative methods.
#' @param maxit Maximum iterations for iterative methods.
#' @param preconditioner `"jacobi"` or `"none"`.
#' @param restart Restart length for GMRES.
#' @param x0 Initial guess for iterative methods.
#' @return `solve_linear_system()` returns a list with `x` and `info`, a list
#'   of diagnostics (`method`, `status`, `residual`, `relative_residual`,
#'   `iterations`, `tol`, `n`, `nnz`, `symmetric`, `time`).
#'   `linear_solver()` returns a list with a function `solve(b)` and `info`.
#' @examples
#' V <- fe_space(pde_mesh_rectangle(nx = 8))
#' A <- assemble_matrix(V, "stiffness") + assemble_matrix(V, "mass")
#' b <- assemble_vector(V, 1)
#' s1 <- solve_linear_system(A, b)
#' s2 <- solve_linear_system(A, b, method = "cg", tol = 1e-12)
#' max(abs(s1$x - s2$x))
#' s2$info$iterations
#' @references Saad, Y. (2003). *Iterative Methods for Sparse Linear
#'   Systems*, 2nd edition. SIAM.
#' @export
solve_linear_system <- function(A, b, method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
                                tol = 1e-10, maxit = 10 * length(b),
                                preconditioner = c("jacobi", "none"), restart = 50,
                                x0 = NULL) {
  method <- match.arg(method)
  preconditioner <- match.arg(preconditioner)
  t0 <- proc.time()[["elapsed"]]
  s <- linear_solver(A, method, tol = tol, maxit = maxit,
                     preconditioner = preconditioner, restart = restart)
  out <- s$solve(b, x0 = x0)
  out$info$time <- proc.time()[["elapsed"]] - t0
  out
}

#' @rdname solve_linear_system
#' @export
linear_solver <- function(A, method = c("auto", "cholesky", "lu", "cg", "bicgstab", "gmres"),
                          tol = 1e-10, maxit = NULL, preconditioner = c("jacobi", "none"),
                          restart = 50) {
  method <- match.arg(method)
  preconditioner <- match.arg(preconditioner)
  A <- methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix")
  n <- nrow(A)
  if (ncol(A) != n) stopf("The system matrix must be square.")
  if (is.null(maxit)) maxit <- 10L * n
  maxit <- check_count(maxit, "maxit", min = 1L)
  restart <- check_count(restart, "restart", min = 1L)
  check_number(tol, "tol", positive = TRUE)
  sym <- Matrix::isSymmetric(A, tol = 1e-12)
  base <- list(n = n, nnz = length(A@x), symmetric = sym, tol = tol)
  factor <- NULL
  note <- NULL
  used <- method
  if (method %in% c("auto", "cholesky")) {
    if (sym) {
      factor <- tryCatch(Matrix::Cholesky(Matrix::forceSymmetric(A, "U"), perm = TRUE, LDL = FALSE),
                         error = function(e) e, warning = function(w) simpleError(conditionMessage(w)))
      if (inherits(factor, "error")) {
        if (method == "cholesky") {
          stopf("Cholesky factorization failed: the matrix is not positive definite.")
        }
        note <- "Cholesky failed (matrix not positive definite); used sparse LU."
        factor <- NULL
        used <- "lu"
      } else {
        used <- "cholesky"
      }
    } else {
      if (method == "cholesky") stopf("Cholesky requires a symmetric matrix.")
      used <- "lu"
    }
  }
  if (used == "lu") {
    factor <- tryCatch(Matrix::lu(A), error = function(e) e)
    if (inherits(factor, "error")) {
      stopf("Sparse LU factorization failed: %s", conditionMessage(factor))
    }
  }
  if (used == "cg" && !sym) stopf("Conjugate gradients require a symmetric matrix.")
  # Pivot ratio of the factorization: a tiny ratio means the matrix is
  # singular to working precision even though the factorization succeeded.
  singular <- FALSE
  if (!is.null(factor)) {
    piv <- tryCatch(abs(as.numeric(if (used == "lu") Matrix::diag(factor@U) else Matrix::diag(factor))),
                    error = function(e) NULL)
    if (length(piv) && max(piv) > 0) {
      ratio <- min(piv) / max(piv)
      singular <- !is.finite(ratio) || ratio < 1e3 * .Machine$double.eps
    }
  }
  dinv <- if (preconditioner == "jacobi") {
    d <- Matrix::diag(A)
    ifelse(d != 0, 1 / d, 1)
  } else rep(1, n)
  solve_fun <- function(b, x0 = NULL) {
    b <- as.numeric(b)
    if (length(b) != n) stopf("The right-hand side has length %d; expected %d.", length(b), n)
    it <- NA_integer_
    if (used %in% c("cholesky", "lu")) {
      x <- tryCatch(as.numeric(Matrix::solve(factor, b)),
                    warning = function(w) stopf("Sparse solve failed: %s", conditionMessage(w)),
                    error = function(e) stopf("Sparse solve failed: %s", conditionMessage(e)))
      status <- "solved"
    } else {
      res <- switch(used,
                    cg = krylov_cg(A, b, x0, dinv, tol, maxit),
                    bicgstab = krylov_bicgstab(A, b, x0, dinv, tol, maxit),
                    gmres = krylov_gmres(A, b, x0, dinv, tol, maxit, restart))
      x <- res$x; it <- res$iterations
      status <- if (res$converged) "converged" else "not converged"
    }
    r <- sqrt(sum((b - as.numeric(A %*% x))^2))
    nb <- sqrt(sum(b^2))
    rel <- if (nb > 0) r / nb else r
    if (used %in% c("cholesky", "lu") && !(is.finite(rel) && rel < 1e-6)) status <- "inaccurate"
    if (singular && status == "solved") status <- "singular"
    if (!all(is.finite(x))) status <- "failed"
    info <- c(list(method = used, status = status, residual = r, relative_residual = rel,
                   iterations = it), base)
    if (!is.null(note)) info$note <- note
    list(x = x, info = info)
  }
  list(solve = solve_fun, info = c(list(method = used), base), matrix = A)
}

krylov_cg <- function(A, b, x0, dinv, tol, maxit) {
  x <- if (is.null(x0)) numeric(length(b)) else as.numeric(x0)
  r <- b - as.numeric(A %*% x)
  nb <- sqrt(sum(b^2)); if (nb == 0) return(list(x = numeric(length(b)), iterations = 0L, converged = TRUE))
  z <- dinv * r; p <- z; rz <- sum(r * z)
  for (k in seq_len(maxit)) {
    if (sqrt(sum(r^2)) <= tol * nb) return(list(x = x, iterations = k - 1L, converged = TRUE))
    Ap <- as.numeric(A %*% p)
    pAp <- sum(p * Ap)
    if (pAp <= 0) return(list(x = x, iterations = k, converged = FALSE))
    a <- rz / pAp
    x <- x + a * p; r <- r - a * Ap
    z <- dinv * r; rz_new <- sum(r * z)
    p <- z + (rz_new / rz) * p; rz <- rz_new
  }
  list(x = x, iterations = maxit, converged = sqrt(sum(r^2)) <= tol * nb)
}

krylov_bicgstab <- function(A, b, x0, dinv, tol, maxit) {
  x <- if (is.null(x0)) numeric(length(b)) else as.numeric(x0)
  r <- b - as.numeric(A %*% x)
  nb <- sqrt(sum(b^2)); if (nb == 0) return(list(x = numeric(length(b)), iterations = 0L, converged = TRUE))
  rhat <- r; rho <- alpha <- omega <- 1
  v <- p <- numeric(length(b))
  for (k in seq_len(maxit)) {
    if (sqrt(sum(r^2)) <= tol * nb) return(list(x = x, iterations = k - 1L, converged = TRUE))
    rho_new <- sum(rhat * r)
    if (rho_new == 0) break
    beta <- (rho_new / rho) * (alpha / omega)
    p <- r + beta * (p - omega * v)
    ph <- dinv * p
    v <- as.numeric(A %*% ph)
    rv <- sum(rhat * v)
    if (!is.finite(rv) || rv == 0) break
    alpha <- rho_new / rv
    s <- r - alpha * v
    if (sqrt(sum(s^2)) <= tol * nb) {
      x <- x + alpha * ph
      return(list(x = x, iterations = k, converged = TRUE))
    }
    sh <- dinv * s
    tt <- as.numeric(A %*% sh)
    t2 <- sum(tt * tt)
    if (!is.finite(t2) || t2 == 0) { x <- x + alpha * ph; r <- s; break }
    omega <- sum(tt * s) / t2
    x <- x + alpha * ph + omega * sh
    r <- s - omega * tt
    rho <- rho_new
    if (omega == 0) break
  }
  list(x = x, iterations = k, converged = isTRUE(sqrt(sum(r^2)) <= tol * nb))
}

krylov_gmres <- function(A, b, x0, dinv, tol, maxit, restart) {
  n <- length(b)
  x <- if (is.null(x0)) numeric(n) else as.numeric(x0)
  nb <- sqrt(sum(b^2)); if (nb == 0) return(list(x = numeric(length(b)), iterations = 0L, converged = TRUE))
  m <- min(restart, n)
  total <- 0L
  while (total < maxit) {
    r <- b - as.numeric(A %*% x)
    beta <- sqrt(sum(r^2))
    if (beta <= tol * nb) return(list(x = x, iterations = total, converged = TRUE))
    V <- matrix(0, n, m + 1L); H <- matrix(0, m + 1L, m)
    cs <- sn <- numeric(m); g <- numeric(m + 1L); g[1] <- beta
    V[, 1] <- r / beta
    j_used <- 0L
    for (j in seq_len(m)) {
      total <- total + 1L
      w <- as.numeric(A %*% (dinv * V[, j]))
      for (i in seq_len(j)) {
        H[i, j] <- sum(w * V[, i]); w <- w - H[i, j] * V[, i]
      }
      H[j + 1L, j] <- sqrt(sum(w^2))
      if (H[j + 1L, j] > 0) V[, j + 1L] <- w / H[j + 1L, j]
      for (i in seq_len(j - 1L)) {
        tmp <- cs[i] * H[i, j] + sn[i] * H[i + 1L, j]
        H[i + 1L, j] <- -sn[i] * H[i, j] + cs[i] * H[i + 1L, j]
        H[i, j] <- tmp
      }
      den <- sqrt(H[j, j]^2 + H[j + 1L, j]^2)
      if (!is.finite(den) || den == 0) break
      cs[j] <- H[j, j] / den; sn[j] <- H[j + 1L, j] / den
      H[j, j] <- den; H[j + 1L, j] <- 0
      g[j + 1L] <- -sn[j] * g[j]; g[j] <- cs[j] * g[j]
      j_used <- j
      if (abs(g[j + 1L]) <= tol * nb || total >= maxit) break
    }
    if (j_used == 0L) break
    y <- backsolve(H[seq_len(j_used), seq_len(j_used), drop = FALSE], g[seq_len(j_used)])
    x <- x + dinv * as.numeric(V[, seq_len(j_used), drop = FALSE] %*% y)
  }
  r <- b - as.numeric(A %*% x)
  list(x = x, iterations = total, converged = isTRUE(sqrt(sum(r^2)) <= tol * nb))
}

# Warn when a linear solve did not reach its tolerance.
warn_solver <- function(info, where = "linear solve") {
  if (!is.null(info$status) && !info$status %in% c("solved", "converged")) {
    warnf("The %s finished with status \"%s\" (relative residual %s).", where, info$status,
          fmt_num(info$relative_residual))
  }
  invisible(info)
}
