# Benchmark protocol. Not part of the package checks (benchmarks/ is in
# .Rbuildignore). Run from the package root:
#
#   Rscript benchmarks/protocol.R [repeats]
#
# Each case runs in a fresh R process started under /usr/bin/time, which
# reports the peak resident set size of that process (`-l` on macOS, `-v`
# on Linux). Three kinds of worker are kept apart:
#
# * "solve": the clean public workload, mesh creation plus one solve_pde()
#   call, and nothing else. Its total time and the peak RSS of its process
#   describe what a user pays for one solve.
# * "components": diagnostic timing of the assembly pieces (local element
#   matrices, sparse construction) on a prepared space, in its own process,
#   so neither its time nor its memory enters the solve figures.
# * "accuracy": one solve followed by the L2 error, in its own process, so
#   the error quadrature enters neither the solve time nor its peak RSS. The
#   error is joined to the solve rows by case and mesh size.
#
# Inside a process the case runs twice to warm up (loading, byte
# compilation, first allocations) and then `repeats` more times serially;
# every repeat is saved, and the median and range are reported. The peak RSS includes
# about 290 MB of R and pkgload::load_all() before any work. Methods are
# compared at a fixed error: for each target, the coarsest mesh of each
# method that reaches it is used.
#
# Output: benchmarks/results/protocol-repeats.csv (one row per repeat),
# benchmarks/results/protocol.csv (summaries) and
# benchmarks/results/protocol-session.txt (commit, R, BLAS, threads, CPU).

args <- commandArgs(trailingOnly = TRUE)

# Worker mode: Rscript protocol.R --worker <kind> <case> <n> <repeats>
if (length(args) && args[1] == "--worker") {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  kind <- args[2]; case <- args[3]; n <- as.integer(args[4]); reps <- as.integer(args[5])
  u <- function(x, y) sin(pi * x) * sin(pi * y)
  f <- function(x, y) 2 * pi^2 * u(x, y)
  degree <- if (grepl("P2", case)) 2 else 1
  method <- if (grepl("amg", case)) "cg" else "auto"
  pc <- if (grepl("amg", case)) amg_preconditioner() else "jacobi"
  el <- function() proc.time()[["elapsed"]]
  if (kind == "solve") {
    run <- function() {
      t0 <- el()
      mesh <- pde_mesh_rectangle(nx = n)
      s <- solve_pde(poisson_problem(mesh, source = f, boundary = dirichlet_bc(0), degree = degree),
                     method = method, preconditioner = pc)
      t1 <- el()
      # Only scalars leave run(), so no solution of an earlier repeat stays
      # alive and inflates the peak RSS.
      out <- list(time = c(total = t1 - t0, assembly = s$diagnostics$assembly_time, solve = s$diagnostics$solve_time),
                  extra = c(error = NA_real_, ndof = s$space$ndof))
      rm(s, mesh)
      out
    }
  } else if (kind == "accuracy") {
    s <- solve_pde(poisson_problem(pde_mesh_rectangle(nx = n), source = f, boundary = dirichlet_bc(0), degree = degree),
                   method = method, preconditioner = pc)
    cat("EXTRA", paste(c("error", "ndof"), c(as.numeric(errornorm(s, u, "L2")), s$space$ndof), sep = "=", collapse = ";"), "\n")
    quit(save = "no")
  } else {
    V <- fe_space(pde_mesh_rectangle(nx = n), degree)
    run <- function() {
      t0 <- el()
      geo <- element_geometry(V, NULL)
      loc <- local_bilinear(geo, diff = coef_qp(1, geo))
      t1 <- el()
      invisible(global_matrix(V, loc))
      t2 <- el()
      list(time = c(total = t2 - t0, local = t1 - t0, sparse = t2 - t1), extra = c(error = NA_real_, ndof = V$ndof))
    }
  }
  invisible(run()); invisible(run())
  r <- lapply(seq_len(reps), function(i) { gc(); run() })
  tm <- do.call(rbind, lapply(r, `[[`, "time"))
  extra <- r[[reps]]$extra
  for (i in seq_len(nrow(tm))) {
    cat("REPEAT", paste(c("repeat", colnames(tm)), c(i, tm[i, ]), sep = "=", collapse = ";"), "\n")
  }
  cat("EXTRA", paste(names(extra), extra, sep = "=", collapse = ";"), "\n")
  quit(save = "no")
}

reps <- if (length(args)) as.integer(args[1]) else 5L
dir.create("benchmarks/results", showWarnings = FALSE)
time_flag <- if (Sys.info()[["sysname"]] == "Darwin") "-l" else "-v"
peak_rss_mb <- function(lines) {
  if (time_flag == "-l") {
    v <- grep("maximum resident set size", lines, value = TRUE)
    as.numeric(sub("^\\s*([0-9]+).*", "\\1", v)) / 2^20         # bytes
  } else {
    v <- grep("Maximum resident set size", lines, value = TRUE)
    as.numeric(sub(".*:\\s*([0-9]+).*", "\\1", v)) / 2^10      # kilobytes
  }
}
parse_kv <- function(line, tag) {
  kv <- strsplit(strsplit(sub(paste0("^", tag, " "), "", trimws(line)), ";")[[1]], "=")
  stats::setNames(suppressWarnings(as.numeric(vapply(kv, `[`, "", 2))), vapply(kv, `[`, "", 1))
}
cases <- list(P1_direct = c(32, 64, 128, 256), P2_direct = c(16, 32, 64, 128), P2_amg = c(16, 32, 64, 128))
reps_rows <- list(); rows <- list()
for (case in names(cases)) for (n in cases[[case]]) for (kind in c("accuracy", "solve", "components")) {
  if (kind == "components" && grepl("amg", case)) next   # same assembly as P2_direct
  if (kind == "accuracy") {
    out <- system2("Rscript", c("benchmarks/protocol.R", "--worker", kind, case, n, 1), stdout = TRUE, stderr = TRUE)
    xl <- grep("^EXTRA", out, value = TRUE)
    if (!length(xl)) { cat(out, sep = "\n"); stop("worker failed: ", kind, " ", case, " ", n) }
    err_now <- parse_kv(xl, "EXTRA")[["error"]]
    next
  }
  out <- suppressWarnings(system2("/usr/bin/time", c(time_flag, "Rscript", "benchmarks/protocol.R", "--worker",
                                                     kind, case, n, reps), stdout = TRUE, stderr = TRUE))
  rl <- grep("^REPEAT", out, value = TRUE); xl <- grep("^EXTRA", out, value = TRUE)
  if (!length(rl) || !length(xl)) { cat(out, sep = "\n"); stop("worker failed: ", kind, " ", case, " ", n) }
  tm <- do.call(rbind, lapply(rl, parse_kv, tag = "REPEAT")); ex <- parse_kv(xl, "EXTRA")
  rss <- peak_rss_mb(out)
  reps_rows[[length(reps_rows) + 1L]] <- data.frame(kind = kind, case = case, n = n, tm, check.names = FALSE)
  med <- apply(tm[, colnames(tm) != "repeat", drop = FALSE], 2, stats::median)
  rows[[length(rows) + 1L]] <- data.frame(kind = kind, case = case, n = n, ndof = ex[["ndof"]], t(med),
                                          total_min = min(tm[, "total"]), total_max = max(tm[, "total"]),
                                          repeats = nrow(tm), error = if (kind == "solve") err_now else NA_real_,
                                          peak_rss_mb = rss)
  cat(sprintf("%-10s %-10s n=%4d ndof=%7d total %.3fs [%.3f, %.3f] rss %.0f MB%s\n", kind, case, n,
              as.integer(ex[["ndof"]]), med[["total"]], min(tm[, "total"]), max(tm[, "total"]), rss,
              if (kind == "solve") sprintf(" err %.2e", err_now) else
                sprintf(" local %.3fs sparse %.3fs", med[["local"]], med[["sparse"]])))
}
fill <- function(l) {
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}
write.csv(fill(reps_rows), "benchmarks/results/protocol-repeats.csv", row.names = FALSE)
res <- fill(rows)
write.csv(res, "benchmarks/results/protocol.csv", row.names = FALSE)

cat("\nClean public solves at fixed L2 error (coarsest mesh of each method reaching the target):\n")
sol <- res[res$kind == "solve", ]
for (target in c(1e-3, 1e-4, 1e-5)) {
  pick <- do.call(rbind, lapply(split(sol, sol$case), function(d) {
    d <- d[d$error <= target, ]
    if (nrow(d)) d[1, c("case", "n", "ndof", "error", "total", "peak_rss_mb")]
  }))
  cat("target", format(target), "\n"); if (!is.null(pick)) print(pick, row.names = FALSE)
}
commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
dirty <- tryCatch(length(system2("git", c("status", "--porcelain", "R"), stdout = TRUE)) > 0, error = function(e) NA)
si <- utils::sessionInfo()
cpu <- if (Sys.info()[["sysname"]] == "Darwin") {
  tryCatch(system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE), error = function(e) NA)
} else tryCatch(sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]), error = function(e) NA)
threads <- Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"))
writeLines(c(paste("commit", commit, if (isTRUE(dirty)) "(uncommitted changes in R/)" else ""),
             R.version.string, paste("Matrix", as.character(packageVersion("Matrix"))),
             paste("BLAS", si$BLAS), paste("LAPACK", si$LAPACK),
             paste("CPU", cpu, "; cores", parallel::detectCores()),
             paste0("threads ", paste(names(threads), ifelse(nzchar(threads), threads, "unset"), sep = "=", collapse = " ")),
             paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
             paste("warmups 2, repeats", reps), format(Sys.time(), tz = "UTC", usetz = TRUE)),
           "benchmarks/results/protocol-session.txt")
