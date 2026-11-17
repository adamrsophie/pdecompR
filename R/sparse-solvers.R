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
#' Iterative methods use Jacobi (diagonal) preconditioning by default;
#' `"ssor"` is symmetric Gauss-Seidel (two sparse triangular solves per
#' application, symmetric for symmetric matrices, so usable with CG); a
#' [multigrid_preconditioner()] gives mesh-independent iteration counts on
#' nested meshes; and any `function(r)` returning an approximation of
#' `A^{-1} r` can be passed. Iterative methods
#' iterate on the scaled system `A y = b / max|b|`, which avoids overflow and
#' underflow. They report `"converged"` only when the true residual of the
#' returned solution satisfies `||b - A x|| <= tol * ||b||`, and
#' `"not converged"` otherwise. Direct methods report the same relative
#' residual, status `"inaccurate"` when the normwise backward error
#' `||b - A x|| / (||A|| ||x|| + ||b||)` exceeds `1e-8`, and status `"near-singular"` when the ratio of the smallest to
#' the largest pivot of the factorization is below about `2e-13`. That ratio
#' is a conditioning diagnostic: the matrix is singular or nearly singular to
#' working precision, for example a stiffness matrix without boundary
#' conditions. A system with no unknowns returns method `"none"`. A small residual says the linear system was solved accurately;
#' it says nothing about the discretization error of the PDE.
#'
#' @param A A square sparse matrix (`Matrix` class) or base matrix.
#' @param b Right-hand side vector, or a matrix of right-hand sides (one per
#'   column; direct methods solve all columns with one factorization).
#' @param method Solver, see Details.
#' @param tol Relative residual tolerance for iterative methods.
#' @param maxit Maximum iterations for iterative methods.
#' @param preconditioner `"jacobi"`, `"ssor"`, `"none"`, a `function(r)`, or
#'   a preconditioner from [multigrid_preconditioner()].
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
                                preconditioner = "jacobi", restart = 50,
                                x0 = NULL) {
  method <- match.arg(method)
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
                          tol = 1e-10, maxit = NULL, preconditioner = "jacobi",
                          restart = 50) {
  method <- match.arg(method)
  if (is.character(preconditioner)) {
    preconditioner <- match.arg(preconditioner, c("jacobi", "ssor", "none"))
  } else if (!is.function(preconditioner) && !inherits(preconditioner, "pde_preconditioner")) {
    stopf("`preconditioner` must be \"jacobi\", \"ssor\", \"none\" or a function.")
  }
  A_res <- if (inherits(A, "Matrix")) A else methods::as(A, "CsparseMatrix")  # for residuals, as passed
  A <- Matrix::drop0(methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix"))
  n <- nrow(A)
  if (ncol(A) != n) stopf("The system matrix must be square.")
  if (n == 0L) {
    # No unknowns (every degree of freedom prescribed): nothing to solve.
    base0 <- list(n = 0L, nnz = 0L, symmetric = TRUE, tol = tol)
    solve0 <- function(b, x0 = NULL) {
      if (length(b) != 0L) stopf("The right-hand side has length %d; expected 0.", length(b))
      list(x = numeric(0), info = c(list(method = "none", status = "solved", residual = 0,
                                         relative_residual = 0, iterations = 0L), base0))
    }
    return(list(solve = solve0, info = c(list(method = "none"), base0), matrix = A))
  }
  if (is.null(maxit)) maxit <- 10L * n
  maxit <- check_count(maxit, "maxit", min = 1L)
  restart <- check_count(restart, "restart", min = 1L)
  check_number(tol, "tol", positive = TRUE)
  sym <- Matrix::isSymmetric(A, tol = 1e-12)
  base <- list(n = n, nnz = length(A@x), symmetric = sym, tol = tol)
  normA <- max(Matrix::rowSums(abs(A)))
  factor <- NULL
  note <- NULL
  used <- method
  # Symmetric diagonal equilibration S A S with S = |diag(A)|^(-1/2) (row
  # norms where the diagonal vanishes) for the direct methods: it preserves
  # symmetry, improves accuracy for badly scaled systems, and makes the
  # pivot-ratio diagnostic independent of diagonal scaling.
  dA <- abs(Matrix::diag(A))
  rn <- Matrix::rowSums(abs(A))
  dA[dA == 0] <- rn[dA == 0]
  dA[dA == 0] <- 1
  sc <- 1 / sqrt(dA)
  Ae <- if (method %in% c("auto", "cholesky", "lu")) Matrix::Diagonal(x = sc) %*% A %*% Matrix::Diagonal(x = sc) else A
  if (method %in% c("auto", "cholesky")) {
    if (sym) {
      factor <- tryCatch(Matrix::Cholesky(Matrix::forceSymmetric(Ae, "U"), perm = TRUE, LDL = FALSE),
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
    factor <- tryCatch(Matrix::lu(Ae), error = function(e) e)
    if (inherits(factor, "error")) {
      stopf("Sparse LU factorization failed: %s", conditionMessage(factor))
    }
  }
  if (used == "cg" && !sym) stopf("Conjugate gradients require a symmetric matrix.")
  # Pivot ratio of the factorization of the equilibrated matrix: a tiny
  # ratio means the matrix is numerically singular or nearly so (a
  # conditioning diagnostic, not a proof of exact singularity) even though
  # the factorization succeeded.
  singular <- FALSE
  if (!is.null(factor)) {
    piv <- tryCatch(abs(as.numeric(if (used == "lu") Matrix::diag(factor@U) else Matrix::diag(factor))),
                    error = function(e) NULL)
    if (length(piv) && max(piv) > 0) {
      ratio <- min(piv) / max(piv)
      singular <- !is.finite(ratio) || ratio < 1e3 * .Machine$double.eps
    }
  }
  pc <- if (!used %in% c("cg", "bicgstab", "gmres")) NULL else make_preconditioner(preconditioner, A)
  solve_fun <- function(b, x0 = NULL) {
    if (is.matrix(b) && ncol(b) != 1L) return(solve_block(b, x0))
    b <- as.numeric(b)
    if (length(b) != n) stopf("The right-hand side has length %d; expected %d.", length(b), n)
    it <- NA_integer_
    if (any(!is.finite(b))) stopf("The right-hand side contains non-finite values.")
    if (used %in% c("cholesky", "lu")) {
      x <- tryCatch(sc * as.numeric(Matrix::solve(factor, sc * b)),
                    warning = function(w) stopf("Sparse solve failed: %s", conditionMessage(w)),
                    error = function(e) stopf("Sparse solve failed: %s", conditionMessage(e)))
      status <- "solved"
    } else {
      # Iterate on b / s, which is exact for a linear system and keeps the
      # recurrences away from overflow and underflow.
      sc <- max(abs(b))
      if (sc == 0) {
        res <- list(x = numeric(n), iterations = 0L, converged = TRUE)
      } else {
        xs <- if (is.null(x0)) NULL else as.numeric(x0) / sc
        res <- switch(used,
                      cg = krylov_cg(A, b / sc, xs, pc, tol, maxit),
                      bicgstab = krylov_bicgstab(A, b / sc, xs, pc, tol, maxit),
                      gmres = krylov_gmres(A, b / sc, xs, pc, tol, maxit, restart))
        res$x <- res$x * sc
      }
      x <- res$x; it <- res$iterations
      status <- if (res$converged) "converged" else "not converged"
    }
    r <- safe_norm(b - as.numeric(A_res %*% x))
    nb <- safe_norm(b)
    rel <- if (nb > 0) r / nb else r
    if (used %in% c("cholesky", "lu")) {
      # Normwise backward error: a backward-stable solve gives rounding level
      # even when ||b|| is tiny compared with ||A|| ||x||.
      den <- normA * safe_norm(x) + nb
      bwe <- if (den > 0) r / den else r
      if (!(is.finite(bwe) && bwe < 1e-8)) status <- "inaccurate"
    }
    # Report convergence only when the true residual meets the tolerance; the
    # recursively updated Krylov residual can drift from it.
    if (status == "converged" && !(is.finite(rel) && rel <= tol * (1 + 1e-6))) status <- "not converged"
    if (singular && status == "solved") status <- "near-singular"
    if (!all(is.finite(x))) status <- "failed"
    info <- c(list(method = used, status = status, residual = r, relative_residual = rel,
                   iterations = it), base)
    if (!is.null(note)) info$note <- note
    list(x = x, info = info)
  }
  # Several right-hand sides: one triangular solve pair for all columns with
  # a factorization, column by column otherwise. The info reports the worst
  # column.
  solve_block <- function(B, X0 = NULL) {
    B <- as.matrix(B)
    if (nrow(B) != n) stopf("The right-hand side has %d rows; expected %d.", nrow(B), n)
    if (any(!is.finite(B))) stopf("The right-hand side contains non-finite values.")
    if (!ncol(B)) return(list(x = matrix(0, n, 0), info = c(list(method = used, status = "solved"), base)))
    if (used %in% c("cholesky", "lu")) {
      X <- tryCatch(sc * as.matrix(Matrix::solve(factor, sc * B)),
                    error = function(e) stopf("Sparse solve failed: %s", conditionMessage(e)))
      R <- B - as.matrix(A_res %*% X)
      r <- apply(R, 2, safe_norm); nb <- apply(B, 2, safe_norm)
      rel <- ifelse(nb > 0, r / nb, r)
      den <- normA * apply(X, 2, safe_norm) + nb
      bwe <- ifelse(den > 0, r / den, r)
      status <- if (!all(is.finite(X))) "failed" else if (!all(is.finite(bwe) & bwe < 1e-8)) "inaccurate"
                else if (singular) "near-singular" else "solved"
      k <- which.max(rel)
      info <- c(list(method = used, status = status, residual = r[k], relative_residual = rel[k],
                     iterations = NA_integer_, columns = ncol(B)), base)
      if (!is.null(note)) info$note <- note
      return(list(x = X, info = info))
    }
    cols <- lapply(seq_len(ncol(B)), function(j) solve_fun(B[, j], if (is.null(X0)) NULL else as.matrix(X0)[, j]))
    rel <- vapply(cols, function(c) c$info$relative_residual, numeric(1))
    ok <- vapply(cols, function(c) solve_ok(c$info), logical(1))
    k <- if (all(ok)) which.max(rel) else which(!ok)[1]
    info <- cols[[k]]$info
    info$iterations <- max(vapply(cols, function(c) c$info$iterations, numeric(1)))
    info$columns <- ncol(B)
    list(x = vapply(cols, function(c) c$x, numeric(n)), info = info)
  }
  list(solve = solve_fun, info = c(list(method = used), base), matrix = A)
}

# Turn the `preconditioner` argument into a function r -> M^{-1} r.
make_preconditioner <- function(pc, A) {
  if (inherits(pc, "pde_preconditioner")) return(pc$setup(A))
  if (is.function(pc)) {
    return(function(r) {
      z <- as.numeric(pc(r))
      if (length(z) != length(r) || any(!is.finite(z))) stopf("The preconditioner returned an invalid vector.")
      z
    })
  }
  switch(pc,
    none = function(r) r,
    jacobi = {
      d <- Matrix::diag(A); dinv <- ifelse(d != 0, 1 / d, 1)
      function(r) dinv * r
    },
    ssor = {
      d <- Matrix::diag(A)
      if (any(d == 0)) stopf("SSOR needs a nonzero diagonal.")
      Lo <- Matrix::tril(A); Up <- Matrix::triu(A)
      function(r) as.numeric(Matrix::solve(Up, d * as.numeric(Matrix::solve(Lo, r))))
    })
}

#' Geometric multigrid preconditioner
#'
#' Builds a V-cycle preconditioner for P1 systems on a hierarchy of nested
#' meshes, each obtained from the previous one by one call of
#' [refine_uniform()] or [refine_mesh()]. Prolongation is exact P1
#' interpolation (new nodes are edge midpoints, taken from the refinement
#' record), coarse operators are Galerkin products `P^T A P`, smoothing is
#' forward Gauss-Seidel before and backward Gauss-Seidel after the coarse
#' correction (so the cycle is symmetric for symmetric `A` and can
#' precondition CG), and the coarsest level is solved directly. With
#' Dirichlet conditions the fine-level prolongation is restricted to the
#' free degrees of freedom automatically.
#'
#' Pass it as `preconditioner` to [solve_linear_system()],
#' [linear_solver()] or [solve_pde()] for elliptic problems, together with
#' an iterative `method`. The iteration count then stays nearly constant as
#' the mesh is refined.
#'
#' @param meshes List of nested meshes, coarse to fine; the last one is the
#'   mesh of the system.
#' @param smoothing Gauss-Seidel sweeps before and after the coarse
#'   correction.
#' @return An object of class `pde_preconditioner`.
#' @examples
#' m0 <- pde_mesh_rectangle(nx = 4)
#' ms <- list(m0)
#' for (k in 1:3) ms[[k + 1]] <- refine_uniform(ms[[k]])
#' prob <- poisson_problem(ms[[4]], source = 1, boundary = dirichlet_bc(0))
#' sol <- solve_pde(prob, method = "cg", preconditioner = multigrid_preconditioner(ms))
#' sol$diagnostics$solver$iterations
#' @export
multigrid_preconditioner <- function(meshes, smoothing = 2L) {
  if (!is.list(meshes) || length(meshes) < 2L || !all(vapply(meshes, inherits, logical(1), "pde_mesh"))) {
    stopf("`meshes` must be a list of at least two nested meshes, coarse to fine.")
  }
  smoothing <- check_count(smoothing, "smoothing", min = 1L)
  P <- lapply(seq_len(length(meshes) - 1L), function(k) refinement_prolongation(meshes[[k]], meshes[[k + 1L]]))
  nf <- nrow(meshes[[length(meshes)]]$nodes)
  structure(list(setup = function(A, free = NULL) {
    n <- nrow(A)
    if (!is.null(free) && max(free) > nf) stopf("The multigrid preconditioner supports P1 spaces only.")
    Pl <- P
    top <- length(Pl)
    if (is.null(free)) {
      if (n != nf) stopf("The system has %d unknowns but the finest mesh has %d nodes; give `free`.", n, nf)
    } else {
      if (length(free) != n) stopf("`free` has %d entries; the system has %d unknowns.", length(free), n)
      # Restrict to the free rows level by level, dropping coarse unknowns
      # that no longer influence any free fine unknown.
      rows <- free
      for (k in top:1) {
        Pt <- Pl[[k]][rows, , drop = FALSE]
        rows <- which(Matrix::colSums(abs(Pt)) > 0)
        Pl[[k]] <- Pt[, rows, drop = FALSE]
      }
    }
    ops <- vector("list", top + 1L)
    ops[[top + 1L]] <- A
    for (k in top:1) ops[[k]] <- Matrix::drop0(Matrix::crossprod(Pl[[k]], ops[[k + 1L]] %*% Pl[[k]]))
    coarse <- linear_solver(ops[[1]])
    Lo <- lapply(ops, Matrix::tril); Up <- lapply(ops, Matrix::triu)
    vcycle <- function(k, b) {
      if (k == 1L) return(coarse$solve(b)$x)
      Ak <- ops[[k]]
      x <- numeric(length(b))
      for (s in seq_len(smoothing)) x <- x + as.numeric(Matrix::solve(Lo[[k]], b - as.numeric(Ak %*% x)))
      r <- b - as.numeric(Ak %*% x)
      x <- x + as.numeric(Pl[[k - 1L]] %*% vcycle(k - 1L, as.numeric(Matrix::crossprod(Pl[[k - 1L]], r))))
      for (s in seq_len(smoothing)) x <- x + as.numeric(Matrix::solve(Up[[k]], b - as.numeric(Ak %*% x)))
      x
    }
    function(r) vcycle(top + 1L, as.numeric(r))
  }, levels = length(meshes)), class = "pde_preconditioner")
}

# P1 prolongation from `coarse` to `fine` = one refinement of `coarse`.
refinement_prolongation <- function(coarse, fine) {
  E <- refinement_edges(coarse, fine); n0 <- nrow(coarse$nodes)
  if (is.null(E)) stopf("Each mesh must come from the previous one by one refine_uniform() or refine_mesh() step.")
  ne <- nrow(E)
  Matrix::sparseMatrix(i = c(seq_len(n0), n0 + rep(seq_len(ne), 2)), j = c(seq_len(n0), E[, 1], E[, 2]),
                       x = c(rep(1, n0), rep(0.5, 2 * ne)), dims = c(n0 + ne, n0))
}

krylov_cg <- function(A, b, x0, pc, tol, maxit) {
  x <- if (is.null(x0)) numeric(length(b)) else as.numeric(x0)
  r <- b - as.numeric(A %*% x)
  nb <- sqrt(sum(b^2)); if (nb == 0) return(list(x = numeric(length(b)), iterations = 0L, converged = TRUE))
  z <- pc(r); p <- z; rz <- sum(r * z)
  for (k in seq_len(maxit)) {
    if (sqrt(sum(r^2)) <= tol * nb) return(list(x = x, iterations = k - 1L, converged = TRUE))
    Ap <- as.numeric(A %*% p)
    pAp <- sum(p * Ap)
    if (pAp <= 0) return(list(x = x, iterations = k, converged = FALSE))
    a <- rz / pAp
    x <- x + a * p; r <- r - a * Ap
    z <- pc(r); rz_new <- sum(r * z)
    p <- z + (rz_new / rz) * p; rz <- rz_new
  }
  list(x = x, iterations = maxit, converged = sqrt(sum(r^2)) <= tol * nb)
}

krylov_bicgstab <- function(A, b, x0, pc, tol, maxit) {
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
    ph <- pc(p)
    v <- as.numeric(A %*% ph)
    rv <- sum(rhat * v)
    if (!is.finite(rv) || rv == 0) break
    alpha <- rho_new / rv
    s <- r - alpha * v
    if (sqrt(sum(s^2)) <= tol * nb) {
      x <- x + alpha * ph
      return(list(x = x, iterations = k, converged = TRUE))
    }
    sh <- pc(s)
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

krylov_gmres <- function(A, b, x0, pc, tol, maxit, restart) {
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
      w <- as.numeric(A %*% pc(V[, j]))
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
    x <- x + pc(as.numeric(V[, seq_len(j_used), drop = FALSE] %*% y))
  }
  r <- b - as.numeric(A %*% x)
  list(x = x, iterations = total, converged = isTRUE(sqrt(sum(r^2)) <= tol * nb))
}

solve_ok <- function(info) !is.null(info$status) && info$status %in% c("solved", "converged")

# Warn when a linear solve did not reach its tolerance.
warn_solver <- function(info, where = "linear solve") {
  if (!is.null(info$status) && !info$status %in% c("solved", "converged")) {
    warnf("The %s finished with status \"%s\" (relative residual %s).", where, info$status,
          fmt_num(info$relative_residual))
  }
  invisible(info)
}

# Euclidean norm without overflow or underflow.
safe_norm <- function(v) {
  m <- max(abs(v))
  if (!is.finite(m) || m == 0) return(m)
  m * sqrt(sum((v / m)^2))
}
