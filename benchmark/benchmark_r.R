# Benchmark R matrix population model analysis
# Uses: popbio, Rage, popdemo packages
# Data: COMPADRE/COMADRE pre-exported as plain RDS
#
# Run: Rscript benchmark/benchmark_r.R

cat("=== R MPM Benchmark ===\n")
cat("Loading packages...\n")

suppressPackageStartupMessages({
  library(popbio)
  library(Rage)
  library(popdemo)
})

# --- Load data ---
load_db <- function(db_name) {
  rds_path <- sprintf("/tmp/%s_plain.rds", db_name)
  if (!file.exists(rds_path)) stop("No data at ", rds_path)
  cat("Loading", db_name, "from", rds_path, "\n")
  d <- readRDS(rds_path)

  n <- length(d$species)
  to_mat <- function(raw, dim) {
    if (is.null(raw) || length(raw) == 0) return(matrix(0, dim, dim))
    actual_dim <- floor(sqrt(length(raw)))
    if (actual_dim^2 == length(raw)) {
      return(matrix(as.numeric(raw), actual_dim, actual_dim))
    }
    matrix(as.numeric(raw), dim, dim)
  }

  A <- vector("list", n)
  U <- vector("list", n)
  Fm <- vector("list", n)
  Cm <- vector("list", n)
  dims <- d$dim

  for (i in 1:n) {
    dim <- ifelse(dims[i] == 0, 1, dims[i])
    A[[i]] <- to_mat(d$matA[[i]], dim)
    U[[i]] <- to_mat(d$matU[[i]], dim)
    Fm[[i]] <- to_mat(d$matF[[i]], dim)
    Cm[[i]] <- to_mat(d$matC[[i]], dim)
    dims[i] <- nrow(A[[i]])
  }

  list(species=d$species, dims=dims, lambda_R=d$lambda,
       A=A, U=U, F=Fm, C=Cm)
}

# --- Benchmarking functions ---

# Lambda: dominant eigenvalue
bench_lambda <- function(A) {
  max(Re(eigen(A)$values))
}

# Sensitivity matrix
bench_sensitivity <- function(A) {
  ev <- eigen(A)
  lv <- eigen(t(A))
  idx <- which.max(Re(ev$values))
  w <- Re(ev$vectors[, idx])
  v <- Re(lv$vectors[, idx])
  w <- w / sum(w)
  v <- v / (sum(v * w))
  outer(v, w)
}

# Elasticity matrix
bench_elasticity <- function(A) {
  S <- bench_sensitivity(A)
  lam <- bench_lambda(A)
  (A * S) / lam
}

# Life expectancy (fundamental matrix approach)
bench_life_expect <- function(U) {
  n <- nrow(U)
  N <- tryCatch(solve(diag(n) - U), error = function(e) NULL)
  if (is.null(N)) return(NA)
  sum(N[1, ])
}

# Net reproductive rate
bench_R0 <- function(U, Fm) {
  n <- nrow(U)
  N <- tryCatch(solve(diag(n) - U), error = function(e) NULL)
  if (is.null(N)) return(NA)
  R0_mat <- Fm %*% N
  R0_mat[1, 1]
}

# Generation time (log(R0)/log(lambda))
bench_gen_time <- function(U, Fm) {
  R0 <- bench_R0(U, Fm)
  if (is.na(R0) || R0 <= 0) return(NA)
  A <- U + Fm
  lam <- bench_lambda(A)
  if (lam <= 0) return(NA)
  log(R0) / log(lam)
}

# Longevity (age at which survivorship drops below threshold)
bench_longevity <- function(U, start = 1, lx_crit = 0.01, xmax = 1000) {
  n <- nrow(U)
  lx <- numeric(xmax + 1)
  lx[1] <- 1.0
  state <- rep(0.0, n)
  state[start] <- 1.0
  for (x in 1:xmax) {
    state <- U %*% state
    lx[x + 1] <- sum(state)
    if (lx[x + 1] < lx_crit) return(x)
  }
  return(xmax)
}

# Vital rate: survival
bench_vr_survival <- function(U) {
  mean(colSums(U))
}

# Full analysis pipeline on a single model
bench_full_analysis <- function(A, U, Fm) {
  lam <- bench_lambda(A)
  sens <- bench_sensitivity(A)
  elas <- bench_elasticity(A)
  le <- bench_life_expect(U)
  R0 <- bench_R0(U, Fm)
  Tg <- bench_gen_time(U, Fm)
  lo <- bench_longevity(U)
  vrs <- bench_vr_survival(U)
  list(lambda=lam, le=le, R0=R0, gen_time=Tg, longevity=lo, vr_surv=vrs)
}

# --- Select representative subset ---
select_models <- function(db, n_sample = 500) {
  valid <- which(!is.na(db$lambda_R) & db$lambda_R > 0.1 & db$lambda_R < 10 &
                 db$dims >= 2 & db$dims <= 20)
  if (length(valid) > n_sample) valid <- valid[1:n_sample]
  valid
}

# --- Time a function over many models ---
time_func <- function(func, ..., n_reps = 3) {
  times <- numeric(n_reps)
  for (r in 1:n_reps) {
    t0 <- proc.time()["elapsed"]
    func(...)
    t1 <- proc.time()["elapsed"]
    times[r] <- t1 - t0
  }
  min(times)
}

# --- Main benchmark ---
run_benchmark <- function(db_name) {
  cat("\n", strrep("=", 60), "\n")
  cat("Benchmarking", toupper(db_name), "\n")
  cat(strrep("=", 60), "\n")

  db <- load_db(db_name)
  ids <- select_models(db)
  cat("Selected", length(ids), "models for benchmarking\n")

  # Warmup
  cat("Warming up...\n")
  for (i in ids[1:min(10, length(ids))]) {
    tryCatch(bench_full_analysis(db$A[[i]], db$U[[i]], db$F[[i]]), error=function(e) NULL)
  }

  # --- Individual function benchmarks ---
  cat("\n--- Individual function benchmarks (all models) ---\n")

  # Lambda
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_lambda(db$A[[i]]), error=function(e) NULL)
  }
  t_lambda <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  lambda:         %.4f s  (%d models, %.1f us/model)\n",
              t_lambda, length(ids), t_lambda / length(ids) * 1e6))

  # Sensitivity
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_sensitivity(db$A[[i]]), error=function(e) NULL)
  }
  t_sens <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  sensitivity:    %.4f s  (%d models, %.1f us/model)\n",
              t_sens, length(ids), t_sens / length(ids) * 1e6))

  # Elasticity
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_elasticity(db$A[[i]]), error=function(e) NULL)
  }
  t_elas <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  elasticity:     %.4f s  (%d models, %.1f us/model)\n",
              t_elas, length(ids), t_elas / length(ids) * 1e6))

  # Life expectancy
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_life_expect(db$U[[i]]), error=function(e) NULL)
  }
  t_le <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  life_expect:    %.4f s  (%d models, %.1f us/model)\n",
              t_le, length(ids), t_le / length(ids) * 1e6))

  # Net reproductive rate
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_R0(db$U[[i]], db$F[[i]]), error=function(e) NULL)
  }
  t_R0 <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  net_repro_rate: %.4f s  (%d models, %.1f us/model)\n",
              t_R0, length(ids), t_R0 / length(ids) * 1e6))

  # Generation time
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_gen_time(db$U[[i]], db$F[[i]]), error=function(e) NULL)
  }
  t_gen <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  gen_time:       %.4f s  (%d models, %.1f us/model)\n",
              t_gen, length(ids), t_gen / length(ids) * 1e6))

  # Longevity
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) tryCatch(bench_longevity(db$U[[i]]), error=function(e) NULL)
  }
  t_long <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  longevity:      %.4f s  (%d models, %.1f us/model)\n",
              t_long, length(ids), t_long / length(ids) * 1e6))

  # Full pipeline
  t0 <- proc.time()["elapsed"]
  for (rep in 1:3) {
    for (i in ids) {
      tryCatch(bench_full_analysis(db$A[[i]], db$U[[i]], db$F[[i]]), error=function(e) NULL)
    }
  }
  t_full <- (proc.time()["elapsed"] - t0) / 3
  cat(sprintf("  full_pipeline:  %.4f s  (%d models, %.1f us/model)\n",
              t_full, length(ids), t_full / length(ids) * 1e6))

  # --- Single matrix benchmarks (amortized over many reps) ---
  cat("\n--- Single matrix benchmarks (amortized over 10000 reps) ---\n")
  sizes <- c(3, 5, 10)
  for (sz in sizes) {
    # Find a model with this size
    idx <- which(db$dims[ids] == sz)
    if (length(idx) == 0) next
    i <- ids[idx[1]]
    A <- db$A[[i]]; U <- db$U[[i]]; Fm <- db$F[[i]]

    n_reps <- 10000

    t0 <- proc.time()["elapsed"]
    for (r in 1:n_reps) bench_lambda(A)
    t_la <- (proc.time()["elapsed"] - t0) / n_reps * 1e6
    cat(sprintf("  lambda %dx%d:      %.1f us/call\n", sz, sz, t_la))

    t0 <- proc.time()["elapsed"]
    for (r in 1:n_reps) bench_sensitivity(A)
    t_se <- (proc.time()["elapsed"] - t0) / n_reps * 1e6
    cat(sprintf("  sensitivity %dx%d: %.1f us/call\n", sz, sz, t_se))

    t0 <- proc.time()["elapsed"]
    for (r in 1:n_reps) bench_full_analysis(A, U, Fm)
    t_fa <- (proc.time()["elapsed"] - t0) / n_reps * 1e6
    cat(sprintf("  full %dx%d:        %.1f us/call\n", sz, sz, t_fa))
  }

  # Return summary
  list(n=length(ids), t_lambda=t_lambda, t_sens=t_sens, t_elas=t_elas,
       t_le=t_le, t_R0=t_R0, t_gen=t_gen, t_long=t_long, t_full=t_full)
}

results <- list()
for (db_name in c("compadre", "comadre")) {
  tryCatch({
    results[[db_name]] <- run_benchmark(db_name)
  }, error = function(e) {
    cat("ERROR (", db_name, "): ", conditionMessage(e), "\n")
  })
}

# --- Save results for comparison ---
cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")
for (db_name in names(results)) {
  r <- results[[db_name]]
  cat(sprintf("%-10s  %d models  full_pipeline=%.4fs (%.1f us/model)\n",
              toupper(db_name), r$n, r$t_full, r$t_full/r$n * 1e6))
}

# Save timing data for comparison
saveRDS(results, "/tmp/benchmark_r_mpm_results.rds")
cat("\nResults saved to /tmp/benchmark_r_mpm_results.rds\n")
