# ckd_msm_M1_cv.R
# ─────────────────────────────────────────────────────────────────────────────
# Repeated 5-fold cross-validation for Model M1 (homogeneous CTMC, CKD progression)
#
# Setup:
#   - 10 repeats × 5 folds = 50 (train, val) splits.
#   - Each repeat regenerates a fresh stratified split by initial CKD stage.
#   - All 50 (repeat, fold) tasks run in parallel via parallel::mclapply.
#   - Per task: fit M1 (with covariates) and M0 (no covariates); record
#     M1 train LL/n, M1 test LL/n, M0 test LL/n; extract 27 parameters from M1;
#     compute M1 expected prevalence on the val set.
#
# Outputs:
#   - M1_cv_metrics.csv : 50 rows + mean/sd summary rows
#   - M1_cv_params.csv  : 50 rows + mean/sd + full_model rows
#   - M1_cv_ll.pdf      : box plot M0 vs M1 (50 points each, mean annotated)
#   - M1_hr_sex.pdf, M1_hr_htn.pdf : forest plots over 50 fits
#   - M1_prevalence.pdf : observed vs predicted prevalences with 50-fit OOS band
#
# Robustness:
#   - compute_prev: strict baseline alignment, LOCF, no look-ahead (unchanged).
#   - Final model loaded from M1_msm_hom.rds (fit on full dataset by ckd_msm_M1.R).

suppressPackageStartupMessages({
  library(msm)
  library(parallel)
})

# ── 0. Setup ──────────────────────────────────────────────────────────────────
cat("CKD M1 — Repeated 5-Fold Cross-Validation (10 × 5 = 50 splits)\n")

script_dir <- tryCatch(
  normalizePath(dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))),
  error = function(e) getwd()
)
data_dir <- file.path(dirname(script_dir), "data_processed")
out_dir  <- Sys.getenv("OUT_DIR", unset = file.path(dirname(script_dir), "results"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out_dir, "M1_cv.log"), split = TRUE)

# Knobs ───────────────────────────────────────────────────────────────────────
N_REPEATS <- as.integer(Sys.getenv("N_REPEATS", unset = "10"))
N_FOLDS   <- 5L
N_TASKS   <- N_REPEATS * N_FOLDS                                 # 50
MC_CORES  <- as.integer(Sys.getenv("MC_CORES",
                                   unset = max(1L, min(parallel::detectCores() - 1L, N_TASKS))))
BASE_SEED <- 42L

cat(sprintf("Tasks: %d  |  parallel workers: %d\n", N_TASKS, MC_CORES))

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
dat <- read.csv(file.path(data_dir, "ckd_panel.csv"), stringsAsFactors = FALSE)
dat$state <- as.integer(dat$state)
dat <- dat[order(dat$patient_id, dat$age), ]
cat(sprintf("  Rows: %d | Patients: %d\n", nrow(dat), length(unique(dat$patient_id))))

# ── 2. Shared model components ────────────────────────────────────────────────
V_init <- rbind(
  c(0,    0.05, 0,    0,    0,    0.005),
  c(0,    0,    0.05, 0,    0,    0.005),
  c(0,    0,    0,    0.05, 0,    0.005),
  c(0,    0,    0,    0,    0.05, 0.005),
  c(0,    0,    0,    0,    0,    0.05 ),
  c(0,    0,    0,    0,    0,    0    )
)

ref_cov   <- list(sex = 0, htn_baseline = 0)
age_times <- seq(40, 85, by = 5)

TRANS <- list(
  list(name="1->2", i=1L, j=2L), list(name="2->3", i=2L, j=3L),
  list(name="3->4", i=3L, j=4L), list(name="4->5", i=4L, j=5L),
  list(name="1->6", i=1L, j=6L), list(name="2->6", i=2L, j=6L),
  list(name="3->6", i=3L, j=6L), list(name="4->6", i=4L, j=6L),
  list(name="5->6", i=5L, j=6L)
)
COV_NAMES <- c("sex", "htn_baseline")

# ── Helper: observed and expected prevalence matrices (unchanged) ─────────────
compute_prev <- function(msm_obj, dat_list, age_times, window = 2.5) {
  n_s     <- 6L
  obs_mat <- matrix(0, nrow = length(age_times), ncol = n_s)
  exp_mat <- matrix(0, nrow = length(age_times), ncol = n_s)

  for (k in seq_along(age_times)) {
    t_target <- age_times[k]

    for (subj in dat_list) {
      past_records <- subj[subj$age < (t_target - window), , drop = FALSE]
      if (nrow(past_records) == 0L) next

      ref <- past_records[nrow(past_records), ]
      t0  <- ref$age
      sk  <- ref$state

      obs_state <- NA
      if (any(subj$state == n_s & subj$age <= t_target)) {
        obs_state <- n_s
      } else {
        if (max(subj$age) < (t_target - window)) next

        near <- subj[abs(subj$age - t_target) <= window & subj$state != n_s, , drop = FALSE]
        if (nrow(near) > 0) {
          obs_state <- near$state[which.min(abs(near$age - t_target))]
        } else {
          obs_state <- sk
        }
      }

      obs_mat[k, obs_state] <- obs_mat[k, obs_state] + 1L

      if (sk == n_s) {
        exp_mat[k, n_s] <- exp_mat[k, n_s] + 1L
      } else {
        dt <- max(t_target - t0, 0)
        P <- tryCatch(
          pmatrix.msm(msm_obj, t = dt, covariates = list(sex=ref$sex, htn_baseline=ref$htn_baseline)),
          error = function(e) NULL
        )
        if (!is.null(P)) exp_mat[k, ] <- exp_mat[k, ] + P[sk, ]
      }
    }

    ot <- sum(obs_mat[k, ]); et <- sum(exp_mat[k, ])
    if (ot > 0) obs_mat[k, ] <- 100 * obs_mat[k, ] / ot
    if (et > 0) exp_mat[k, ] <- 100 * exp_mat[k, ] / et
  }
  rownames(obs_mat) <- rownames(exp_mat) <- as.character(age_times)
  colnames(obs_mat) <- colnames(exp_mat) <- paste0("State ", 1:n_s)
  list(obs = obs_mat, exp = exp_mat)
}

# ── Helper: held-out log-likelihood per transition ────────────────────────────
compute_llpern <- function(msm_obj, dat_list) {
  ll_total <- 0
  n_trans  <- 0L

  for (subj in dat_list) {
    if (nrow(subj) < 2L) next
    cov_l <- list(sex = subj$sex[1], htn_baseline = subj$htn_baseline[1])
    for (r in seq_len(nrow(subj) - 1L)) {
      dt <- subj$age[r + 1L] - subj$age[r]
      if (dt <= 0) next
      P <- tryCatch(
        pmatrix.msm(msm_obj, t = dt, covariates = cov_l),
        error = function(e) NULL
      )
      if (is.null(P)) next
      p_ij <- P[subj$state[r], subj$state[r + 1L]]
      ll_total <- ll_total + log(max(p_ij, 1e-300))
      n_trans  <- n_trans + 1L
    }
  }
  if (n_trans == 0L) return(NA_real_)
  ll_total / n_trans
}

# ── Helper: extract parameters from a fitted msm object ───────────────────────
extract_params <- function(msm_obj, ref_cov, TRANS, COV_NAMES) {
  V_ref   <- qmatrix.msm(msm_obj, covariates = ref_cov)$estimates
  hr_list <- hazard.msm(msm_obj)

  params <- c()
  for (tr in TRANS) {
    params[paste0("a_", tr$name)] <- V_ref[tr$i, tr$j]

    trans_key <- sprintf("State %d - State %d", tr$i, tr$j)
    for (cv in COV_NAMES) {
      if (cv %in% names(hr_list)) {
        hr_mat <- as.matrix(hr_list[[cv]])
        if (trans_key %in% rownames(hr_mat)) {
          hr_val <- hr_mat[trans_key, "HR"]
          hr_lo  <- hr_mat[trans_key, "L"]
          hr_hi  <- hr_mat[trans_key, "U"]
        } else {
          hr_val <- hr_lo <- hr_hi <- NA_real_
        }
      } else {
        hr_val <- hr_lo <- hr_hi <- NA_real_
      }
      params[paste0("HR_",         tr$name, "_", cv)] <- hr_val
      params[paste0("beta_",       tr$name, "_", cv)] <- log(hr_val)
      params[paste0("beta_lower_", tr$name, "_", cv)] <- log(hr_lo)
      params[paste0("beta_upper_", tr$name, "_", cv)] <- log(hr_hi)
    }
  }
  params
}

# ── 3. Generate 10 stratified 5-fold partitions ───────────────────────────────
cat("\nGenerating 10 stratified 5-fold partitions...\n")
pt_df <- dat[!duplicated(dat$patient_id), c("patient_id", "sex", "htn_baseline", "state")]
names(pt_df)[names(pt_df) == "state"] <- "init_state"

# fold_assign[[r]] is a named integer vector: patient_id -> fold (1..5)
fold_assign <- vector("list", N_REPEATS)
for (r in seq_len(N_REPEATS)) {
  set.seed(BASE_SEED + r)
  assignment <- integer(nrow(pt_df))
  for (s in unique(pt_df$init_state)) {
    idx <- which(pt_df$init_state == s)
    idx <- sample(idx)
    assignment[idx] <- ((seq_along(idx) - 1L) %% N_FOLDS) + 1L
  }
  names(assignment)  <- as.character(pt_df$patient_id)
  fold_assign[[r]]   <- assignment
}

# Quick sanity check on the first partition
cat("Repeat 1 fold × initial-stage distribution:\n")
print(table(Fold = fold_assign[[1]], InitStage = pt_df$init_state))

# Group dat by patient_id once
full_dat_list <- split(dat[order(dat$patient_id, dat$age), ], dat$patient_id)

# ── 4. Worker: fit M1 + M0 for one (repeat, fold) pair ────────────────────────
fit_one <- function(task) {
  r <- task$repeat_idx
  k <- task$fold

  assignment <- fold_assign[[r]]
  pids       <- names(assignment)

  train_ids <- pids[assignment != k]
  val_ids   <- pids[assignment == k]

  train_dat <- dat[as.character(dat$patient_id) %in% train_ids, ]
  train_dat <- train_dat[order(train_dat$patient_id, train_dat$age), ]
  val_sub   <- dat[as.character(dat$patient_id) %in% val_ids, ]
  val_sub   <- val_sub[order(val_sub$patient_id, val_sub$age), ]

  train_dat_list <- split(train_dat, train_dat$patient_id)
  val_dat_list   <- split(val_sub,   val_sub$patient_id)

  result <- list(
    repeat_idx      = r,
    fold            = k,
    n_train         = length(train_ids),
    n_val           = length(val_ids),
    llpern_m1_train = NA_real_,
    llpern_m1       = NA_real_,
    llpern_m0       = NA_real_,
    params          = NULL,
    exp_val         = NULL,
    m1_ok           = FALSE,
    m0_ok           = FALSE,
    msg             = ""
  )

  # Try a sequence of inits, stop on first success. Failed folds usually have
  # data that makes msm's likelihood either overflow (init too large) or
  # underflow (init too small) — sweep a few scales of V_init to find a safe one.
  fit_msm <- function(use_cov) {
    cov_arg <- if (use_cov) ~ sex + htn_baseline else NULL
    ctl     <- list(fnscale = 4000, maxit = 1500, reltol = 1e-7)

    attempts <- list(
      list(label = "gen.inits=TRUE",   qmat = V_init,       gen = TRUE),
      list(label = "V_init x1",        qmat = V_init,       gen = FALSE),
      list(label = "V_init x2",        qmat = V_init * 2,   gen = FALSE),
      list(label = "V_init x0.5",      qmat = V_init * 0.5, gen = FALSE),
      list(label = "V_init x5",        qmat = V_init * 5,   gen = FALSE),
      list(label = "V_init x0.2",      qmat = V_init * 0.2, gen = FALSE)
    )

    msgs <- character(0)
    for (att in attempts) {
      res <- tryCatch(
        msm(state ~ age, subject = patient_id, data = train_dat,
            qmatrix = att$qmat, covariates = cov_arg, gen.inits = att$gen,
            control = ctl),
        error = function(e) e
      )
      if (!inherits(res, "error")) {
        msg <- if (length(msgs) > 0L)
          sprintf("OK on '%s' after: %s", att$label, paste(msgs, collapse = " | "))
        else ""
        return(list(fit = res, msg = msg))
      }
      msgs <- c(msgs, sprintf("%s: %s", att$label, conditionMessage(res)))
    }
    list(fit = NULL, msg = paste(msgs, collapse = " | "))
  }

  m1_attempt <- fit_msm(use_cov = TRUE)
  msm_m1 <- m1_attempt$fit
  if (nchar(m1_attempt$msg) > 0L) result$msg <- paste0("M1: ", m1_attempt$msg)

  if (!is.null(msm_m1)) {
    result$m1_ok           <- TRUE
    result$llpern_m1_train <- compute_llpern(msm_m1, train_dat_list)
    result$llpern_m1       <- compute_llpern(msm_m1, val_dat_list)
    result$exp_val         <- compute_prev(msm_m1, val_dat_list, age_times)$exp
    result$params          <- tryCatch(
      extract_params(msm_m1, ref_cov, TRANS, COV_NAMES),
      error = function(e) NULL
    )
  }

  m0_attempt <- fit_msm(use_cov = FALSE)
  msm_m0 <- m0_attempt$fit
  if (nchar(m0_attempt$msg) > 0L)
    result$msg <- paste0(result$msg, if (nchar(result$msg) > 0L) " || " else "", "M0: ", m0_attempt$msg)

  if (!is.null(msm_m0)) {
    result$m0_ok     <- TRUE
    result$llpern_m0 <- compute_llpern(msm_m0, val_dat_list)
  }

  result
}

# ── 5. Build task list and run in parallel ────────────────────────────────────
tasks <- vector("list", N_TASKS)
ix <- 1L
for (r in seq_len(N_REPEATS)) {
  for (k in seq_len(N_FOLDS)) {
    tasks[[ix]] <- list(repeat_idx = r, fold = k)
    ix <- ix + 1L
  }
}

cat(sprintf("\nRunning %d tasks across %d workers...\n", N_TASKS, MC_CORES))
t_start <- Sys.time()
RNGkind("L'Ecuyer-CMRG")
set.seed(BASE_SEED)

results <- mclapply(tasks, fit_one, mc.cores = MC_CORES, mc.preschedule = FALSE)

t_end <- Sys.time()
cat(sprintf("Parallel fitting done in %.1f s (wall clock).\n",
            as.numeric(difftime(t_end, t_start, units = "secs"))))

# Brief summary of any failures
n_m1_ok <- sum(sapply(results, function(x) isTRUE(x$m1_ok)))
n_m0_ok <- sum(sapply(results, function(x) isTRUE(x$m0_ok)))
cat(sprintf("M1 converged: %d / %d   |   M0 converged: %d / %d\n",
            n_m1_ok, N_TASKS, n_m0_ok, N_TASKS))
for (res in results) {
  if (nchar(res$msg) > 0L)
    cat(sprintf("  (repeat=%d, fold=%d) %s\n", res$repeat_idx, res$fold, res$msg))
}

# ── 6. Aggregate metrics and parameters ───────────────────────────────────────
cv_metrics_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    repeat_idx      = x$repeat_idx,
    fold            = x$fold,
    llpern_m1_train = x$llpern_m1_train,
    llpern_m1       = x$llpern_m1,
    llpern_m0       = x$llpern_m0,
    stringsAsFactors = FALSE
  )
}))
cv_metrics_df <- cv_metrics_df[order(cv_metrics_df$repeat_idx, cv_metrics_df$fold), ]
rownames(cv_metrics_df) <- NULL

cat("\nPer-task metrics (first 10 rows):\n")
print(head(cv_metrics_df, 10))

cv_params_list <- lapply(results, function(x) {
  if (is.null(x$params)) return(NULL)
  c(repeat_idx = x$repeat_idx, fold = x$fold, x$params)
})
cv_params_list <- Filter(Negate(is.null), cv_params_list)
cv_params_df <- if (length(cv_params_list) > 0L) {
  do.call(rbind, lapply(cv_params_list,
                        function(x) as.data.frame(t(x), check.names = FALSE,
                                                  stringsAsFactors = FALSE)))
} else NULL
if (!is.null(cv_params_df)) {
  cv_params_df <- cv_params_df[order(cv_params_df$repeat_idx, cv_params_df$fold), ]
  rownames(cv_params_df) <- NULL
}

# Expected-prevalence array: (n_age_times × 6 × n_valid)
exp_val_list <- lapply(results, function(x) x$exp_val)
exp_val_list <- Filter(Negate(is.null), exp_val_list)

# ── 7. Load final model (fit on full dataset by ckd_msm_M1.R) ─────────────────
rds_path <- Sys.getenv("MSM_FINAL_RDS", unset = file.path(out_dir, "M1_msm_hom.rds"))
if (!file.exists(rds_path)) stop("M1_msm_hom.rds not found in results/. Run ckd_msm_M1.R first.")
cat(sprintf("\nLoading final model from %s\n", rds_path))
msm_final <- readRDS(rds_path)

final_params  <- extract_params(msm_final, ref_cov, TRANS, COV_NAMES)
final_hr_list <- hazard.msm(msm_final)
full_llpern   <- compute_llpern(msm_final, full_dat_list)

cat("Computing full-data observed and expected prevalences...\n")
prev_full      <- compute_prev(msm_final, full_dat_list, age_times)
obs_full       <- prev_full$obs
exp_full_final <- prev_full$exp

# ── 8. Save CSV outputs ───────────────────────────────────────────────────────
metric_num_cols <- setdiff(names(cv_metrics_df), c("repeat_idx", "fold"))
mean_m <- sapply(metric_num_cols, function(c) round(mean(cv_metrics_df[[c]], na.rm = TRUE), 6))
sd_m   <- sapply(metric_num_cols, function(c) round(sd  (cv_metrics_df[[c]], na.rm = TRUE), 6))

metrics_out <- cv_metrics_df
metrics_out$repeat_idx <- as.character(metrics_out$repeat_idx)
metrics_out$fold       <- as.character(metrics_out$fold)
mean_row <- as.data.frame(c(list(repeat_idx = "mean", fold = ""), as.list(mean_m)), stringsAsFactors = FALSE)
sd_row   <- as.data.frame(c(list(repeat_idx = "sd",   fold = ""), as.list(sd_m)),   stringsAsFactors = FALSE)
metrics_out <- rbind(metrics_out, mean_row, sd_row)
write.csv(metrics_out, file.path(out_dir, "M1_cv_metrics.csv"), row.names = FALSE)
cat(sprintf("Saved: %s  (%d rows + mean/sd)\n",
            file.path(out_dir, "M1_cv_metrics.csv"), nrow(cv_metrics_df)))

if (!is.null(cv_params_df)) {
  param_num_cols <- setdiff(names(cv_params_df), c("repeat_idx", "fold"))
  mean_p <- sapply(param_num_cols, function(p) round(mean(cv_params_df[[p]], na.rm = TRUE), 6))
  sd_p   <- sapply(param_num_cols, function(p) round(sd  (cv_params_df[[p]], na.rm = TRUE), 6))
  full_p <- round(final_params, 6)

  params_out <- cv_params_df
  params_out$repeat_idx <- as.character(params_out$repeat_idx)
  params_out$fold       <- as.character(params_out$fold)

  make_summary_row <- function(label, vec) {
    row <- c(list(repeat_idx = label, fold = ""), as.list(vec))
    as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  }
  mean_pr <- make_summary_row("mean",       mean_p)
  sd_pr   <- make_summary_row("sd",         sd_p)
  full_pr <- make_summary_row("full_model", full_p)

  # Align column order/names explicitly to avoid rbind name-mismatch errors
  mean_pr <- mean_pr[, names(params_out), drop = FALSE]
  sd_pr   <- sd_pr  [, names(params_out), drop = FALSE]
  full_pr <- full_pr[, names(params_out), drop = FALSE]
  params_out <- rbind(params_out, mean_pr, sd_pr, full_pr)

  write.csv(params_out, file.path(out_dir, "M1_cv_params.csv"), row.names = FALSE)
  cat(sprintf("Saved: %s  (%d rows + mean/sd + full)\n",
              file.path(out_dir, "M1_cv_params.csv"), nrow(cv_params_df)))
}

# ── 9. Plots ──────────────────────────────────────────────────────────────────
state_colors <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b")
state_names  <- c("Stage 1","Stage 2","Stage 3","Stage 4","Stage 5","Death")

# ── Plot 1: Box plot — M0 vs M1 test LL/n across 50 splits ────────────────────
{
  ll_plot_path <- file.path(out_dir, "M1_cv_ll.pdf")
  pdf(ll_plot_path, width = 5.1, height = 5.2, pointsize = 9)
  par(mar = c(4.8, 4.5, 1.2, 1.2), mgp = c(2.5, 0.7, 0))

  m0_vals <- cv_metrics_df$llpern_m0
  m1_vals <- cv_metrics_df$llpern_m1
  m0_finite <- m0_vals[is.finite(m0_vals)]
  m1_finite <- m1_vals[is.finite(m1_vals)]

  m0_mean <- mean(m0_finite)
  m1_mean <- mean(m1_finite)
  m0_sd   <- sd(m0_finite)
  m1_sd   <- sd(m1_finite)

  y_all  <- c(m0_finite, m1_finite, full_llpern)
  y_pad  <- diff(range(y_all)) * 0.15
  y_lim  <- c(min(y_all) - y_pad, max(y_all) + y_pad)

  box_col <- c("#d62728", "#2ca02c")  # M0 red, M1 green

  n_splits <- length(m1_finite)
  x_labels <- c(
    sprintf("M0 (null)\n%.4f ± %.4f",        m0_mean, m0_sd),
    sprintf("M1 (sex + HTN)\n%.4f ± %.4f",   m1_mean, m1_sd)
  )

  bp <- boxplot(list(M0 = m0_finite, M1 = m1_finite),
                outline = FALSE, ylim = y_lim,
                col = adjustcolor(box_col, alpha.f = 0.18),
                border = box_col, boxwex = 0.55, whisklty = 1, staplewex = 0.4,
                ylab = "Test log-likelihood per transition (LL/n)",
                xaxt = "n")
  axis(1, at = c(1, 2), labels = x_labels, padj = 0.6, cex.axis = 0.9, tick = FALSE)

  set.seed(BASE_SEED)
  points(jitter(rep(1, length(m0_finite)), amount = 0.12), m0_finite,
         pch = 21, bg = adjustcolor(box_col[1], alpha.f = 0.55),
         col = box_col[1], cex = 0.85)
  points(jitter(rep(2, length(m1_finite)), amount = 0.12), m1_finite,
         pch = 21, bg = adjustcolor(box_col[2], alpha.f = 0.55),
         col = box_col[2], cex = 0.85)

  points(c(1, 2), c(m0_mean, m1_mean), pch = 23, bg = "white",
         col = box_col, cex = 1.6, lwd = 1.8)

  abline(h = full_llpern, lty = 2, col = "grey35", lwd = 1.2)
  text(1.5, full_llpern,
       sprintf("Full-data M1\n%.4f", full_llpern),
       pos = 3, offset = 0.3, cex = 0.8, col = "grey25")

  legend("bottomright",
         legend = c(sprintf("%d test LL/n (%d × %d-fold)",
                            n_splits, N_REPEATS, N_FOLDS),
                    "Mean"),
         pch = c(21, 23),
         pt.bg = c(adjustcolor("grey60", 0.55), "white"),
         col   = c("grey40", "black"),
         pt.cex = c(0.9, 1.4),
         bty = "n", cex = 0.82)

  dev.off()
  cat(sprintf("Plot saved to: %s\n", ll_plot_path))

  cat(sprintf("  M0 test LL/n: mean = %.5f, sd = %.5f (n = %d)\n",
              m0_mean, m0_sd, length(m0_finite)))
  cat(sprintf("  M1 test LL/n: mean = %.5f, sd = %.5f (n = %d)\n",
              m1_mean, m1_sd, length(m1_finite)))
}

# ── Plot 2: Forest plots — HR on log-scaled axis (over 50 fits) ───────────────
if (!is.null(cv_params_df)) {
  trans_order <- sapply(TRANS, `[[`, "name")

  for (cv_covar in COV_NAMES) {
    hr_fname     <- if (cv_covar == "sex") "M1_hr_sex.pdf" else "M1_hr_htn.pdf"
    hr_plot_path <- file.path(out_dir, hr_fname)
    hr_title     <- if (cv_covar == "sex") "(a) Sex (male vs. female)"
                    else                   "(b) Hypertension (HTN vs. no HTN)"
    pdf(hr_plot_path, width = 3.5, height = 4.5, pointsize = 9)
    par(mar = c(4.0, 6.5, 1.8, 1.2), mgp = c(2.5, 0.7, 0))

    hr_pnames <- paste0("HR_", trans_order, "_", cv_covar)

    full_hr <- unlist(final_params[hr_pnames])
    full_lo <- full_hi <- numeric(length(TRANS))
    for (idx in seq_along(TRANS)) {
      tr     <- TRANS[[idx]]
      msmkey <- sprintf("State %d - State %d", tr$i, tr$j)
      if (msmkey %in% names(final_hr_list)) {
        hr_mat <- as.matrix(final_hr_list[[msmkey]])
        if (cv_covar %in% rownames(hr_mat)) {
          full_lo[idx] <- hr_mat[cv_covar, "L"]
          full_hi[idx] <- hr_mat[cv_covar, "U"]
        }
      }
    }

    cv_log_means <- sapply(hr_pnames, function(p) mean(log(cv_params_df[[p]]), na.rm = TRUE))
    cv_log_sds   <- sapply(hr_pnames, function(p) sd  (log(cv_params_df[[p]]), na.rm = TRUE))
    cv_geo_means <- exp(cv_log_means)
    cv_geo_lo    <- exp(cv_log_means - cv_log_sds)
    cv_geo_hi    <- exp(cv_log_means + cv_log_sds)

    n_tr    <- length(TRANS)
    row_lbl <- trans_order
    x_lim   <- c(0.3, 35)

    plot(NULL, xlim = x_lim, ylim = c(0.5, n_tr + 0.5),
         log  = "x",
         xlab = "Hazard ratio (log scale)", ylab = "",
         main = hr_title, font.main = 1, cex.main = 1,
         yaxt = "n")
    axis(2, at = seq_len(n_tr), labels = rev(row_lbl), las = 1, cex.axis = 0.85)
    abline(v = 1, lty = 2, col = "grey55", lwd = 1.2)

    for (idx in seq_len(n_tr)) {
      y_hi <- n_tr - idx + 1 + 0.13
      y_lo <- n_tr - idx + 1 - 0.13

      segments(full_lo[idx], y_hi, full_hi[idx], y_hi, col = "black", lwd = 1.5)
      points(full_hr[idx], y_hi, pch = 19, col = "black", cex = 1.0)

      segments(cv_geo_lo[idx], y_lo, cv_geo_hi[idx], y_lo, col = "#1f77b4", lwd = 1.5)
      points(cv_geo_means[idx], y_lo, pch = 23, col = "#1f77b4", bg = "white", cex = 1.0)
    }

    legend("topright",
           legend = c("Full model", "10×5-fold CV mean ± 1 SD"),
           pch = c(19, 23), col = c("black", "#1f77b4"), pt.bg = c(NA, "white"),
           lty = c(NA, 1), lwd = 1.5, bty = "n", cex = 0.75)
    dev.off()
    cat(sprintf("Plot saved to: %s\n", hr_plot_path))
  }
}

# ── Plot 3: Observed vs Expected prevalence (OOS = average of 50 val-set fits)
{
  n_valid <- length(exp_val_list)

  if (n_valid > 0) {
    prev_plot_path <- file.path(out_dir, "M1_prevalence.pdf")
    pdf(prev_plot_path, width = 5.1, height = 5.0, pointsize = 9)

    fold_arr <- array(unlist(exp_val_list), dim = c(length(age_times), 6L, n_valid))
    cv_mean  <- apply(fold_arr, c(1, 2), mean, na.rm = TRUE)
    cv_sd    <- apply(fold_arr, c(1, 2), sd,   na.rm = TRUE)

    col_full <- "#1f77b4"
    col_oos  <- "#ff7f0e"
    par(mfrow = c(2, 3), mar = c(4.0, 3.8, 2.5, 1.2), oma = c(0, 0, 0.5, 0), mgp = c(2.2, 0.6, 0))

    for (s in 1:6) {
      obs_v <- obs_full[, s]
      exp_v <- exp_full_final[, s]
      cv_v  <- cv_mean[, s]

      cv_lo <- pmax(cv_v - cv_sd[, s], 0)
      cv_hi <- cv_v + cv_sd[, s]

      y_max <- max(c(obs_v, exp_v, cv_hi), na.rm = TRUE) * 1.15
      if (!is.finite(y_max) || y_max == 0) y_max <- 1
      if (s == 1) y_max <- y_max * 2.2

      plot(NULL, xlim = range(age_times), ylim = c(0, y_max),
           xlab = "Age (years)", ylab = "% patients", main = state_names[s],
           font.main = 1, cex.main = 1)

      polygon(c(age_times, rev(age_times)), c(cv_hi, rev(cv_lo)),
              col = adjustcolor(col_oos, alpha.f = 0.20), border = NA)

      lines(age_times, exp_v, lty = 1, col = col_full, lwd = 2.0)
      lines(age_times, cv_v,  lty = 2, col = col_oos,  lwd = 1.5)
      points(age_times, obs_v, pch = 16, col = "black", cex = 0.8)

      if (s == 1) {
        legend("topright",
               legend = c("Observed", "Predicted (100%)", "Predicted (OOS, 10×5)", "± 1 SD band"),
               lty  = c(NA, 1, 2, NA),
               lwd  = c(NA, 2.0, 1.5, NA),
               pch  = c(16, NA, NA, 15),
               col  = c("black", col_full, col_oos, adjustcolor(col_oos, alpha.f = 0.4)),
               bty  = "n", cex = 0.80, pt.cex = c(0.8, 1, 1, 1.5))
      }
    }
    dev.off()
    cat(sprintf("Plot saved to: %s\n", prev_plot_path))
  }
}

cat("===== Done =====\n")
sink()
