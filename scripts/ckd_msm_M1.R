# ckd_msm_M1.R
# ─────────────────────────────────────────────────────────────────────────────
# M1: Strictly homogeneous CTMC for CKD progression
#
# Model:  lambda_ij(z_i) = a_ij * exp(beta_ij^T z_i)
#         z_i is FIXED per patient (no time-varying covariates)
#         => V does not depend on time t => strict homogeneous CTMC
#
# States:
#   1 = CKD Stage 1 (eGFR >= 90)
#   2 = CKD Stage 2 (eGFR 60-89)
#   3 = CKD Stage 3 (eGFR 30-59)
#   4 = CKD Stage 4 (eGFR 15-29)
#   5 = CKD Stage 5 (eGFR < 15)
#   6 = Death (absorbing)
#
# Allowed transitions: 1->2, 2->3, 3->4, 4->5, 1->6, 2->6, 3->6, 4->6, 5->6
# No backward transitions, no stage-skipping.
#
# Covariates (all time-invariant, fixed at first observation per patient):
#   sex          : 1 = male, 0 = female
#   htn_baseline : 1 = hypertension at CKD diagnosis, 0 = no
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages(library(msm))

# ── 0. Sink setup (must be before any cat/print) ──────────────────────────────
script_dir <- tryCatch(
  normalizePath(dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))),
  error = function(e) getwd()
)
data_dir <- file.path(dirname(script_dir), "data_processed")
out_dir  <- file.path(dirname(script_dir), "results")
dir.create(out_dir, showWarnings = FALSE)
sink(file.path(out_dir, "M1.log"), split = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
dat <- read.csv(file.path(data_dir, "ckd_panel.csv"), stringsAsFactors = FALSE)

cat(sprintf("  Rows: %d | Patients: %d\n", nrow(dat), length(unique(dat$patient_id))))
cat("  State counts:\n"); print(table(dat$state))

# ── 2. Preprocessing ──────────────────────────────────────────────────────────
dat$state <- as.integer(dat$state)
dat       <- dat[order(dat$patient_id, dat$age), ]

# ── Observed transition counts ────────────────────────────────────────────────
cat("\n--- Observed Transition Counts ---\n")
dat_tc            <- dat
dat_tc$next_state <- c(dat_tc$state[-1], NA_integer_)
dat_tc$next_id    <- c(dat_tc$patient_id[-1], NA_character_)
valid_tc          <- dat_tc[!is.na(dat_tc$next_id) & dat_tc$patient_id == dat_tc$next_id, ]
trans_tbl         <- table(From = valid_tc$state, To = valid_tc$next_state)
cat("Transition count matrix (rows=From, cols=To):\n")
print(trans_tbl)
allowed_arcs <- list(c(1,2),c(2,3),c(3,4),c(4,5),c(1,6),c(2,6),c(3,6),c(4,6),c(5,6))
tc_df <- data.frame(
  transition = sapply(allowed_arcs, function(x) sprintf("%d->%d", x[1], x[2])),
  count      = sapply(allowed_arcs, function(x) {
    v <- trans_tbl[as.character(x[1]), as.character(x[2])]
    if (length(v) > 0 && !is.na(v)) as.integer(v) else 0L
  }),
  stringsAsFactors = FALSE
)
cat("\nAllowed M1 arcs — observed counts:\n")
print(tc_df, row.names = FALSE)

# Per-transition breakdown by sex and hypertension
cat("\n--- Transition breakdown by sex and hypertension ---\n")
trans_detail <- do.call(rbind, lapply(allowed_arcs, function(arc) {
  sub <- valid_tc[valid_tc$state == arc[1] & valid_tc$next_state == arc[2], ]
  data.frame(
    transition = sprintf("%d->%d", arc[1], arc[2]),
    total      = nrow(sub),
    male       = sum(sub$sex == 1),
    female     = sum(sub$sex == 0),
    htn_yes    = sum(sub$htn_baseline == 1),
    htn_no     = sum(sub$htn_baseline == 0),
    stringsAsFactors = FALSE
  )
}))
print(trans_detail, row.names = FALSE)

# Cohort summary
n_pts   <- length(unique(dat$patient_id))
pt_1row <- dat[!duplicated(dat$patient_id), ]   
n_male  <- sum(pt_1row$sex == 1)
n_fem   <- sum(pt_1row$sex == 0)
n_htn   <- sum(pt_1row$htn_baseline == 1)
n_htn_m <- sum(pt_1row$sex == 1 & pt_1row$htn_baseline == 1)
n_htn_f <- sum(pt_1row$sex == 0 & pt_1row$htn_baseline == 1)

cat("\n--- Cohort Summary ---\n")
cat(sprintf("  Total patients : %d\n", n_pts))
cat(sprintf("  Male           : %d (%.1f%%)\n", n_male, 100 * n_male / n_pts))
cat(sprintf("  Female         : %d (%.1f%%)\n", n_fem,  100 * n_fem  / n_pts))
cat(sprintf("  HTN at CKD dx  : %d (%.1f%%)  — Male: %d, Female: %d\n",
            n_htn, 100 * n_htn / n_pts, n_htn_m, n_htn_f))

# ── 3. Initial V matrix ───────────────────────────────────────────────────────
V_init <- rbind(
  c(0,    0.05, 0,    0,    0,    0.005),  
  c(0,    0,    0.05, 0,    0,    0.005),  
  c(0,    0,    0,    0.05, 0,    0.005),  
  c(0,    0,    0,    0,    0.05, 0.005),  
  c(0,    0,    0,    0,    0,    0.05 ),  
  c(0,    0,    0,    0,    0,    0    )   
)

# ── 4. Fit M1 ─────────────────────────────────────────────────────────────────
cat("\nFitting M1...\n")
msm_hom <- msm(
  state ~ age,
  subject    = patient_id,
  data       = dat,
  qmatrix    = V_init,
  covariates = ~ sex + htn_baseline,
  gen.inits  = TRUE,  
  control    = list(fnscale = 4000, maxit = 3000, reltol = 1e-8)
)
cat("M1 converged.\n")

# ── 5. Full parameter table ───────────────────────────────────────────────────
ref_cov <- list(sex = 0, htn_baseline = 0)

cat("\n\n====== M1: Complete Parameter Table ======\n")
cat("Model: lambda_ij(z) = a_ij * exp(beta_ij^T z)\n")
cat("Reference group: female, no HTN\n")

# (a) Baseline rates
cat("\n--- Baseline rates a_ij (at reference group) ---\n")
V_ref <- qmatrix.msm(msm_hom, covariates = ref_cov)
n <- nrow(V_ref$estimates)
base_rows <- data.frame()
for (i in seq_len(n)) for (j in seq_len(n)) {
  if (i != j && V_ref$estimates[i, j] > 1e-10) {
    base_rows <- rbind(base_rows, data.frame(
      transition = sprintf("%d->%d", i, j),
      a_ij       = round(V_ref$estimates[i, j], 6),
      CI_lower   = round(V_ref$L[i, j],         6),
      CI_upper   = round(V_ref$U[i, j],         6)
    ))
  }
}
print(base_rows, row.names = FALSE)

# (b) Covariate coefficients
cat("\n--- Covariate effects: beta_ij = log(HR_ij), HR_ij = exp(beta_ij) ---\n\n")
hr_list  <- hazard.msm(msm_hom)
cov_rows <- data.frame()
for (cv_covar in names(hr_list)) {
  hr_mat <- as.matrix(hr_list[[cv_covar]])
  for (trans in rownames(hr_mat)) {
    hr <- hr_mat[trans, "HR"]
    lo <- hr_mat[trans, "L"]
    hi <- hr_mat[trans, "U"]
    cov_rows <- rbind(cov_rows, data.frame(
      transition = trans,
      covariate  = cv_covar,
      beta       = round(log(hr), 4),
      beta_lower = round(log(lo), 4),
      beta_upper = round(log(hi), 4),
      HR         = round(hr, 4),
      HR_lower   = round(lo, 4),
      HR_upper   = round(hi, 4)
    ))
  }
}
print(cov_rows, row.names = FALSE)

# ── 6. Mean sojourn times: 4 profiles ────────────────────────────────────────
cat("\n--- Mean Sojourn Times (4 profiles: sex x HTN) ---\n\n")
covprofiles <- expand.grid(sex = c(0, 1), htn_baseline = c(0, 1))
covprofiles$label <- with(covprofiles, sprintf("sex=%d htn=%d", sex, htn_baseline))

sojourn_table <- data.frame()
for (i in seq_len(nrow(covprofiles))) {
  cov_list <- list(sex = covprofiles$sex[i], htn_baseline = covprofiles$htn_baseline[i])
  lbl     <- covprofiles$label[i]
  soj     <- sojourn.msm(msm_hom, covariates = cov_list)
  soj_est <- soj$estimates
  soj_lo  <- soj$L
  soj_hi  <- soj$U
  row_df <- data.frame(profile = lbl, sex = covprofiles$sex[i], htn_baseline = covprofiles$htn_baseline[i])
  for (k in seq_along(soj_est)) {
    sname <- paste0("State", k)
    row_df[[paste0(sname, "_est")]] <- round(soj_est[k], 3)
    row_df[[paste0(sname, "_lo")]]  <- round(soj_lo[k],  3)
    row_df[[paste0(sname, "_hi")]]  <- round(soj_hi[k],  3)
  }
  sojourn_table <- rbind(sojourn_table, row_df)
}
print(sojourn_table)

# ── 7. Validation (Strict Formula 5 Alignment) ────────────────────────────────
cat("\n\n========== VALIDATION ==========\n")

cat("\n--- Pearson goodness-of-fit (M1) ---\n")
tryCatch(print(pearson.msm(msm_hom, timegroups = 4)),
         error = function(e) cat("Pearson test failed:", conditionMessage(e), "\n"))

age_times <- seq(40, 85, by = 5)
cat("\n--- Observed vs Expected state prevalences (Formula 5 method) ---\n")

compute_prev_formula5 <- function(msm_obj, dat, age_times, window = 2.5) {
  n_s     <- 6L
  obs_mat <- matrix(0, nrow = length(age_times), ncol = n_s)
  exp_mat <- matrix(0, nrow = length(age_times), ncol = n_s)
  
  dat      <- dat[order(dat$patient_id, dat$age), ]
  dat_list <- split(dat, dat$patient_id)

  for (k in seq_along(age_times)) {
    t_target <- age_times[k]
    
    for (subj in dat_list) {
      # 1. 寻找预测起点 t0 (必须在目标时间窗口之前)
      past_records <- subj[subj$age < (t_target - window), , drop = FALSE]
      if (nrow(past_records) == 0L) next 
      
      ref <- past_records[nrow(past_records), ]
      t0  <- ref$age
      sk  <- ref$state
      
      # 2. 确定 t_target 时的观测状态
      obs_state <- NA
      
      # 情形 A：在 t_target 之前已经确认死亡 (严格不穿越未来)
      if (any(subj$state == n_s & subj$age <= t_target)) {
        obs_state <- n_s
      } else {
        # 情形 B：检查是否已经失访 (最后一次出现都在窗口前)
        if (max(subj$age) < (t_target - window)) {
          next  # 剔除失联者，对齐分母
        }
        
        # 情形 C：寻找窗口内的观测 (严防未来死亡记录污染)
        near <- subj[abs(subj$age - t_target) <= window & subj$state != n_s, , drop = FALSE]
        if (nrow(near) > 0) {
          obs_state <- near$state[which.min(abs(near$age - t_target))]
        } else {
          # 情形 D：窗口内没来，结转最后一次已知状态
          obs_state <- sk
        }
      }
      
      # 3. 累加人数
      obs_mat[k, obs_state] <- obs_mat[k, obs_state] + 1L
      
      if (sk == n_s) {
        exp_mat[k, n_s] <- exp_mat[k, n_s] + 1L
      } else {
        dt <- max(t_target - t0, 0) 
        P <- tryCatch(
          pmatrix.msm(msm_obj, t = dt,
                      covariates = list(sex          = ref$sex,
                                        htn_baseline = ref$htn_baseline)),
          error = function(e) NULL
        )
        if (!is.null(P)) {
          exp_mat[k, ] <- exp_mat[k, ] + P[sk, ]
        }
      }
    }
    
    # 4. 转为百分比
    ot <- sum(obs_mat[k, ])
    et <- sum(exp_mat[k, ])
    if (ot > 0) obs_mat[k, ] <- 100 * obs_mat[k, ] / ot
    if (et > 0) exp_mat[k, ] <- 100 * exp_mat[k, ] / et
  }
  
  rownames(obs_mat) <- rownames(exp_mat) <- as.character(age_times)
  colnames(obs_mat) <- colnames(exp_mat) <- paste0("State ", 1:n_s)
  list(`Observed percentages` = obs_mat, `Expected percentages` = exp_mat)
}

prev_hom <- compute_prev_formula5(msm_hom, dat, age_times)
cat("Observed percentages:\n");  print(round(prev_hom$`Observed percentages`,  2))
cat("Expected percentages:\n");  print(round(prev_hom$`Expected percentages`,  2))

# ── 8. Plots ──────────────────────────────────────────────────────────────────
profile_cols   <- c("#1f77b4","#1f77b4","#ff7f0e","#ff7f0e")
profile_ltys   <- c(1, 2, 1, 2)
profile_labels <- c("Female, no HTN", "Male, no HTN", "Female, HTN", "Male, HTN")
surv_mat <- matrix(NA, nrow = nrow(covprofiles), ncol = 5,
                   dimnames = list(covprofiles$label, paste0("Stage ", 1:5)))
for (i in seq_len(nrow(covprofiles))) {
  Q    <- qmatrix.msm(msm_hom,
                      covariates = list(sex          = covprofiles$sex[i],
                                        htn_baseline = covprofiles$htn_baseline[i]))$estimates
  Q_TT <- Q[1:5, 1:5]
  N    <- solve(-Q_TT)
  surv_mat[i, ] <- rowSums(N)
}

# Plot 1: Mean survival time from each starting state
pdf(file.path(out_dir, "M1_survival_time.pdf"), width = 3.1, height = 3.5, pointsize = 9)
par(mar = c(4.0, 4.0, 1.8, 1.2), mgp = c(2.2, 0.6, 0))
plot(NULL, xlim = c(1, 5), ylim = c(0, max(surv_mat) * 1.1),
     xlab = "Starting CKD stage", ylab = "Expected years until death",
     main = "(a) Mean survival time", font.main = 1, cex.main = 1, xaxt = "n")
axis(1, at = 1:5, labels = paste0("Stage ", 1:5))
for (i in seq_len(nrow(covprofiles))) {
  lines(1:5, surv_mat[i, ], col = profile_cols[i], lty = profile_ltys[i],
        lwd = 1.5, type = "b", pch = 16, cex = 0.8)
}
legend("bottomleft", legend = profile_labels,
       col = profile_cols, lty = profile_ltys, lwd = 1.5, pch = 16,
       pt.cex = 0.8, bty = "n", cex = 1)
dev.off()

# Plot 2: P(alive | start Stage 1)
t_seq <- seq(0, 30, by = 0.5)
pdf(file.path(out_dir, "M1_survival_curve.pdf"), width = 3.1, height = 3.5, pointsize = 9)
par(mar = c(4.0, 4.0, 1.8, 1.2), mgp = c(2.2, 0.6, 0))
plot(NULL, xlim = c(0, 30), ylim = c(0, 1),
     xlab = "Years since CKD Stage 1 diagnosis",
     ylab = "Probability of being alive",
     main = "(b) Survival probability from Stage 1", font.main = 1, cex.main = 1)
abline(h = 0.5, lty = 3, col = "grey70")
for (i in seq_len(nrow(covprofiles))) {
  p_alive <- sapply(t_seq, function(tt)
    1 - pmatrix.msm(msm_hom, t = tt,
          covariates = list(sex          = covprofiles$sex[i],
                            htn_baseline = covprofiles$htn_baseline[i]))[1, 6])
  lines(t_seq, p_alive, lwd = 1.5,
        col = profile_cols[i], lty = profile_ltys[i])
}
legend("bottomleft", legend = profile_labels,
       col = profile_cols, lty = profile_ltys, lwd = 1.5, bty = "n", cex = 1)
dev.off()

cat(sprintf("Plots saved to: %s\n", file.path(out_dir, "M1_survival_*.pdf")))

# ── 9. Save results ───────────────────────────────────────────────────────────
write.csv(round(V_ref$estimates, 6), file.path(out_dir, "M1_V_reference.csv"))
write.csv(round(qmatrix.msm(msm_hom)$estimates, 6), file.path(out_dir, "M1_V_mean.csv"))
write.csv(base_rows,     file.path(out_dir, "M1_params_baseline.csv"),  row.names = FALSE)
write.csv(cov_rows,      file.path(out_dir, "M1_params_covariates.csv"), row.names = FALSE)
write.csv(sojourn_table, file.path(out_dir, "M1_sojourn_4profiles.csv"), row.names = FALSE)
write.csv(round(prev_hom$`Observed percentages`, 2), file.path(out_dir, "M1_prevalence_observed.csv"))
write.csv(round(prev_hom$`Expected percentages`, 2), file.path(out_dir, "M1_prevalence_expected.csv"))
saveRDS(msm_hom, file.path(out_dir, "M1_msm_hom.rds"))

cat(sprintf("\nAll results saved to: %s\n", normalizePath(out_dir)))
sink()