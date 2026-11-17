# Benchmark protocol. Not part of the package checks (benchmarks/ is in
# .Rbuildignore). Run from the package root:
#
#   Rscript benchmarks/protocol.R [repeats]
#
# Every case runs in a fresh R process started under /usr/bin/time, which
# reports the true peak resident set size of that process (`-l` on macOS,
# `-v` on Linux). Inside the process the case runs once to warm up (loading,
# byte compilation, first allocations) and then `repeats` more times
# serially; the median and the range of those repeats are reported.
# Methods are compared at a fixed error: for each target, the coarsest mesh
# of each method that reaches it is used. Results go to
# benchmarks/results/protocol.csv with the session information.

args <- commandArgs(trailingOnly = TRUE)

# Worker mode: Rscript protocol.R --worker <case> <n> <repeats>
if (length(args) && args[1] == "--worker") {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  case <- args[2]; n <- as.integer(args[3]); reps <- as.integer(args[4])
  u <- function(x, y) sin(pi * x) * sin(pi * y)
  f <- function(x, y) 2 * pi^2 * u(x, y)
  degree <- if (grepl("P2", case)) 2 else 1
  method <- if (grepl("amg", case)) "cg" else "auto"
  pc <- if (grepl("amg", case)) amg_preconditioner() else "jacobi"
  run <- function() {
    t0 <- proc.time()[["elapsed"]]
    mesh <- pde_mesh_rectangle(nx = n)
    V <- fe_space(mesh, degree)
    ta <- proc.time()[["elapsed"]]
    geo <- element_geometry(V, NULL)
    loc <- local_bilinear(geo, diff = coef_qp(1, geo))
    tl <- proc.time()[["elapsed"]]
    invisible(global_matrix(V, loc))
    tg <- proc.time()[["elapsed"]]
    s <- solve_pde(poisson_problem(mesh, source = f, boundary = dirichlet_bc(0), degree = degree),
                   method = method, preconditioner = pc)
    t1 <- proc.time()[["elapsed"]]
    c(total = t1 - t0, local = tl - ta, sparse = tg - tl, assembly = s$diagnostics$assembly_time,
      solve = s$diagnostics$solve_time, error = as.numeric(errornorm(s, u, "L2")), ndof = s$space$ndof)
  }
  invisible(run())
  r <- sapply(seq_len(reps), function(i) run())
  out <- c(apply(r[c("total", "local", "sparse", "assembly", "solve"), , drop = FALSE], 1, stats::median),
           total_min = min(r["total", ]), total_max = max(r["total", ]), error = unname(r["error", 1]), ndof = unname(r["ndof", 1]))
  cat("RESULT", paste(names(out), out, sep = "=", collapse = ";"), "\n")
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
cases <- list(P1_direct = c(32, 64, 128, 256), P2_direct = c(16, 32, 64, 128), P2_amg = c(16, 32, 64, 128))
rows <- list()
for (case in names(cases)) for (n in cases[[case]]) {
  out <- suppressWarnings(system2("/usr/bin/time", c(time_flag, "Rscript", "benchmarks/protocol.R", "--worker",
                                                     case, n, reps), stdout = TRUE, stderr = TRUE))
  res <- grep("^RESULT", out, value = TRUE)
  if (!length(res)) { cat(out, sep = "\n"); stop("worker failed: ", case, " ", n) }
  kv <- strsplit(strsplit(sub("^RESULT ", "", trimws(res)), ";")[[1]], "=")
  val <- stats::setNames(as.numeric(vapply(kv, `[`, "", 2)), vapply(kv, `[`, "", 1))
  rows[[length(rows) + 1L]] <- data.frame(case = case, n = n, t(val), peak_rss_mb = peak_rss_mb(out))
  cat(sprintf("%-10s n=%4d ndof=%7d total %.3fs [%.3f, %.3f] local %.3fs sparseMatrix %.3fs solve %.3fs rss %.0f MB err %.2e\n",
              case, n, val[["ndof"]], val[["total"]], val[["total_min"]], val[["total_max"]], val[["local"]],
              val[["sparse"]], val[["solve"]], rows[[length(rows)]]$peak_rss_mb, val[["error"]]))
}
res <- do.call(rbind, rows)
res$sparse_share <- res$sparse / (res$local + res$sparse)
write.csv(res, "benchmarks/results/protocol.csv", row.names = FALSE)

cat("\nAt fixed L2 error (coarsest mesh of each method reaching the target):\n")
for (target in c(1e-3, 1e-4, 1e-5)) {
  pick <- do.call(rbind, lapply(split(res, res$case), function(d) {
    d <- d[d$error <= target, ]
    if (nrow(d)) d[1, c("case", "n", "ndof", "error", "total", "peak_rss_mb")]
  }))
  cat("target", format(target), "\n"); if (!is.null(pick)) print(pick, row.names = FALSE)
}
writeLines(c(R.version.string, paste("Matrix", as.character(packageVersion("Matrix"))),
             paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
             paste("repeats", reps), format(Sys.time(), tz = "UTC", usetz = TRUE)),
           "benchmarks/results/protocol-session.txt")
