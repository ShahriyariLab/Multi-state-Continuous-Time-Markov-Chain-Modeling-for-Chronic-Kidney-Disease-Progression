# ckd_msm_M1_sensitivity.R
# ─────────────────────────────────────────────────────────────────────────────
# Sensitivity analysis for M1: htn_baseline beta(1->6) constrained to CI bounds
# to assess how this wide-CI parameter propagates uncertainty to the rest of
# the model. Since the 95% CI excludes 0, only plausible extremes are tested.
#
#   S_lb: htn_baseline beta(1->6) fixed at lower 95% CI bound of full model
#   S_ub: htn_baseline beta(1->6) fixed at upper 95% CI bound of full model
#
# CI bounds are read automatically from results/M1_params_covariates.csv.
#
# Requires: results/M1_msm_hom.rds           (from ckd_msm_M1.R)
#           results/M1_params_covariates.csv   (from ckd_msm_M1.R)
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(msm))

cat("CKD M1 — Sensitivity: htn_baseline beta(1->6) constrained (S_lb / S_ub)\n")

# ── 0. Paths ──────────────────────────────────────────────────────────────────
script_dir <- tryCatch(
  normalizePath(dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))),
  error = function(e) getwd()
)
data_dir <- file.path(dirname(script_dir), "data_processed")
out_dir  <- file.path(dirname(script_dir), "results")
dir.create(out_dir, showWarnings = FALSE)
sink(file.path(out_dir, "M1_sensitivity.log"), split = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
cat("Loading data...\n")
dat <- read.csv(file.path(data_dir, "ckd_panel.csv"), stringsAsFactors = FALSE)
dat$state <- as.integer(dat$state)
dat       <- dat[order(dat$patient_id, dat$age), ]
cat(sprintf("  Rows: %d | Patients: %d\n", nrow(dat), length(unique(dat$patient_id))))

# ── 2. Shared components ──────────────────────────────────────────────────────
V_init <- rbind(
  c(0,    0.05, 0,    0,    0,    0.005),
  c(0,    0,    0.05, 0,    0,    0.005),
  c(0,    0,    0,    0.05, 0,    0.005),
  c(0,    0,    0,    0,    0.05, 0.005),
  c(0,    0,    0,    0,    0,    0.05 ),
  c(0,    0,    0,    0,    0,    0    )
)
ref_cov   <- list(sex = 0, htn_baseline = 0)
COV_NAMES <- c("sex", "htn_baseline")
TRANS_Q   <- list(c(1,2),c(1,6),c(2,3),c(2,6),c(3,4),c(3,6),c(4,5),c(4,6),c(5,6))
TRANS_NM  <- c("1->2","1->6","2->3","2->6","3->4","3->6","4->5","4->6","5->6")

# Full: both covariates on every transition
COV_FULL <- list(
  "1-2" = ~ sex + htn_baseline, "1-6" = ~ sex + htn_baseline,
  "2-3" = ~ sex + htn_baseline, "2-6" = ~ sex + htn_baseline,
  "3-4" = ~ sex + htn_baseline, "3-6" = ~ sex + htn_baseline,
  "4-5" = ~ sex + htn_baseline, "4-6" = ~ sex + htn_baseline,
  "5-6" = ~ sex + htn_baseline
)
# ── Helpers ───────────────────────────────────────────────────────────────────
extract_all <- function(msm_obj, ref_cov, TRANS_Q, TRANS_NM, COV_NAMES) {
  V_ref   <- qmatrix.msm(msm_obj, covariates = ref_cov)
  hr_list <- hazard.msm(msm_obj)

  rate_rows <- data.frame(
    block = "baseline_rate", transition = TRANS_NM, covariate = NA_character_,
    a_ij    = sapply(TRANS_Q, function(x) V_ref$estimates[x[1], x[2]]),
    a_lower = sapply(TRANS_Q, function(x) V_ref$L[x[1], x[2]]),
    a_upper = sapply(TRANS_Q, function(x) V_ref$U[x[1], x[2]]),
    beta = NA_real_, beta_lower = NA_real_, beta_upper = NA_real_,
    HR   = NA_real_, HR_lower   = NA_real_, HR_upper   = NA_real_,
    significant = NA, fixed = FALSE,
    stringsAsFactors = FALSE
  )

  cov_rows <- do.call(rbind, lapply(COV_NAMES, function(cv) {
    hr_mat <- if (cv %in% names(hr_list)) as.matrix(hr_list[[cv]]) else NULL
    do.call(rbind, lapply(seq_along(TRANS_Q), function(k) {
      trans_key <- sprintf("State %d - State %d", TRANS_Q[[k]][1], TRANS_Q[[k]][2])
      if (!is.null(hr_mat) && trans_key %in% rownames(hr_mat)) {
        hr  <- hr_mat[trans_key, "HR"]
        lo  <- hr_mat[trans_key, "L"]
        hi  <- hr_mat[trans_key, "U"]
        fix <- FALSE
        sig <- !is.na(lo) && !(lo <= 1 & hi >= 1)
      } else {
        hr <- 1; lo <- 1; hi <- 1; fix <- TRUE; sig <- FALSE
      }
      data.frame(
        block = "covariate", transition = TRANS_NM[k], covariate = cv,
        a_ij = NA_real_, a_lower = NA_real_, a_upper = NA_real_,
        beta = log(hr), beta_lower = if (fix) 0 else log(lo),
        beta_upper = if (fix) 0 else log(hi),
        HR = hr, HR_lower = lo, HR_upper = hi,
        significant = sig, fixed = fix,
        stringsAsFactors = FALSE
      )
    }))
  }))

  rbind(rate_rows, cov_rows)
}

print_params <- function(df, label) {
  cat(sprintf("\n%s\n  %s\n%s\n", strrep("─", 90), label, strrep("─", 90)))
  r <- df[df$block == "baseline_rate", ]
  cat(sprintf("  %-8s  %12s  (%10s, %10s)\n", "trans", "a_ij", "lower", "upper"))
  for (i in seq_len(nrow(r)))
    cat(sprintf("  %-8s  %12.6f  (%10.6f, %10.6f)\n",
                r$transition[i], r$a_ij[i], r$a_lower[i], r$a_upper[i]))
  cat(sprintf("\n  %-8s  %-14s  %8s  %18s  %8s  %18s  %s\n",
              "trans","covariate","beta","beta 95% CI","HR","HR 95% CI","status"))
  cv <- df[df$block == "covariate", ]
  for (i in seq_len(nrow(cv))) {
    tag <- if (cv$fixed[i]) "[FIXED]" else if (cv$significant[i]) "[*]" else "[ ]"
    cat(sprintf("  %-8s  %-14s  %8.4f  (%7.4f,%7.4f)  %8.4f  (%7.4f,%7.4f)  %s\n",
                cv$transition[i], cv$covariate[i],
                cv$beta[i], cv$beta_lower[i], cv$beta_upper[i],
                cv$HR[i],   cv$HR_lower[i],   cv$HR_upper[i], tag))
  }
}

compare_cov <- function(df_full, df_cstr, sfx) {
  f <- df_full[df_full$block=="covariate",
               c("transition","covariate","beta","HR","HR_lower","HR_upper","significant","fixed")]
  c <- df_cstr[df_cstr$block=="covariate",
               c("transition","covariate","beta","HR","HR_lower","HR_upper","significant","fixed")]
  names(f)[3:8] <- paste0(names(f)[3:8], "_full")
  names(c)[3:8] <- paste0(names(c)[3:8], "_", sfx)
  m <- merge(f, c, by = c("transition","covariate"))
  m$delta_beta  <- round(m[[paste0("beta_",sfx)]] - m$beta_full, 4)
  m$pct_change  <- round(100 * m$delta_beta / pmax(abs(m$beta_full), 1e-6), 2)
  m$sig_changed <- m$significant_full != m[[paste0("significant_",sfx)]] &
                   !m[[paste0("fixed_",sfx)]]
  m
}

print_comparison <- function(m, sfx) {
  cat(sprintf("\n  %-8s %-14s %9s %9s %8s %7s  %5s %7s  %s\n",
              "trans","covariate","beta_full",paste0("beta_",sfx),
              "delta","pct%","sig_f",paste0("sig_",sfx),"note"))
  for (i in seq_len(nrow(m))) {
    tag <- if (isTRUE(m[[paste0("fixed_",sfx)]][i])) "(fixed)" else
           if (isTRUE(m$sig_changed[i])) "<< SIGNIFICANCE CHANGED" else ""
    cat(sprintf("  %-8s %-14s %9.4f %9.4f %8.4f %6.1f%%  %5s %7s  %s\n",
                m$transition[i], m$covariate[i],
                m$beta_full[i], m[[paste0("beta_",sfx)]][i],
                m$delta_beta[i], m$pct_change[i],
                m$significant_full[i], m[[paste0("significant_",sfx)]][i], tag))
  }
}

do_lrt <- function(msm_constrained, msm_ref, df) {
  chi2 <- msm_constrained$minus2loglik - msm_ref$minus2loglik
  p    <- pchisq(chi2, df = df, lower.tail = FALSE)
  list(ll_ref  = -msm_ref$minus2loglik / 2,
       ll_cstr = -msm_constrained$minus2loglik / 2,
       chi2 = chi2, df = df, p = p)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Load full model + read CI bounds for htn_baseline(1->6)
# ═══════════════════════════════════════════════════════════════════════════════
rds_path     <- file.path(out_dir, "M1_msm_hom.rds")
cov_csv_path <- file.path(out_dir, "M1_params_covariates.csv")
if (!file.exists(rds_path))     stop("M1_msm_hom.rds not found. Run ckd_msm_M1.R first.")
if (!file.exists(cov_csv_path)) stop("M1_params_covariates.csv not found. Run ckd_msm_M1.R first.")

cat("\nLoading full model from M1_msm_hom.rds...\n")
msm_full    <- readRDS(rds_path)
params_full <- extract_all(msm_full, ref_cov, TRANS_Q, TRANS_NM, COV_NAMES)
print_params(params_full, "Full Model (M1)")

params_cov <- read.csv(cov_csv_path, stringsAsFactors = FALSE)
htn16_row  <- params_cov[params_cov$transition == "State 1 - State 6" & params_cov$covariate == "htn_baseline", ]
beta_full  <- htn16_row$beta
beta_lb    <- htn16_row$beta_lower
beta_ub    <- htn16_row$beta_upper

cat(sprintf("\nhtn_baseline(1->6) from full model:\n"))
cat(sprintf("  beta = %.4f  (95%% CI: %.4f, %.4f)\n", beta_full, beta_lb, beta_ub))
cat(sprintf("  HR   = %.4f  (95%% CI: %.4f, %.4f)\n", exp(beta_full), exp(beta_lb), exp(beta_ub)))

# ── Parameter index + covinits setup for fixedpars ───────────────────────────
# msm() does NOT have an `inits` argument. The correct way to pass initial
# values for covariate betas is via `covinits` (named list, one vector per
# covariate, length = number of allowed transitions).
# fixedpars then fixes that parameter at the value supplied in covinits.
#
# Parameter vector order in msm opt$par:
#   [1..9]   log baseline rates  (one per allowed transition)
#   [10..18] sex betas           (one per transition, same order)
#   [19..27] htn_baseline betas  (one per transition, same order)
#   1->6 is the 2nd transition → htn_baseline(1->6) = index 9+9+2 = 20
full_pars  <- msm_full$opt$par
par_names  <- names(full_pars)
n_trans    <- length(TRANS_Q)                      # 9
htn16_idx  <- grep("htn_baseline.*1-6", par_names)
if (length(htn16_idx) == 0) {
  htn16_idx <- n_trans + n_trans + 2L              # fallback: 9+9+2 = 20
  cat(sprintf("  Parameter index (computed): %d\n", htn16_idx))
} else {
  cat(sprintf("  Parameter index (from names): %d — '%s'\n", htn16_idx, par_names[htn16_idx]))
}
cat(sprintf("  Current value in opt$par[%d]: %.6f\n", htn16_idx, full_pars[htn16_idx]))

# Extract per-transition beta vectors from full model to use as warm start
sex_betas  <- full_pars[(n_trans + 1L):(2L * n_trans)]   # indices 10-18
htn_betas  <- full_pars[(2L * n_trans + 1L):(3L * n_trans)]  # indices 19-27
htn16_pos  <- 2L   # 1->6 is 2nd in TRANS_Q

# Use full model's baseline rates as starting Q matrix (better convergence)
V_full_ref <- qmatrix.msm(msm_full, covariates = ref_cov)$estimates
cat(sprintf("  sex betas (full):          %s\n", paste(round(sex_betas, 4), collapse=", ")))
cat(sprintf("  htn_baseline betas (full): %s\n", paste(round(htn_betas, 4), collapse=", ")))

# ═══════════════════════════════════════════════════════════════════════════════
# S_lb: htn_baseline beta(1->6) fixed at lower 95% CI bound
# ═══════════════════════════════════════════════════════════════════════════════
cat(sprintf("\n\n====== S_lb: htn_baseline beta(1->6) fixed at lower CI = %.4f (HR=%.4f) ======\n",
            beta_lb, exp(beta_lb)))

# covinits: same as full model for all betas, except htn_baseline(1->6) = beta_lb
htn_covinits_lb <- htn_betas
htn_covinits_lb[htn16_pos] <- beta_lb

msm_slb <- tryCatch(
  msm(state ~ age, subject = patient_id, data = dat,
      qmatrix   = V_full_ref, covariates = COV_FULL, gen.inits = FALSE,
      covinits  = list(sex = sex_betas, htn_baseline = htn_covinits_lb),
      fixedpars = htn16_idx,
      control   = list(fnscale = 4000, maxit = 3000, reltol = 1e-8)),
  error = function(e) { cat("ERROR S_lb:", conditionMessage(e), "\n"); NULL }
)

params_slb <- NULL
if (!is.null(msm_slb)) {
  cat("S_lb converged.\n")
  params_slb <- extract_all(msm_slb, ref_cov, TRANS_Q, TRANS_NM, COV_NAMES)
  print_params(params_slb, sprintf("S_lb: htn_baseline beta(1->6) fixed = %.4f (lower CI)", beta_lb))
  cmp_slb <- compare_cov(params_full, params_slb, "Slb")
  cat("\nParameter changes (full -> S_lb):\n"); print_comparison(cmp_slb, "Slb")
  write.csv(cmp_slb,    file.path(out_dir, "M1_Slb_htn16lb_comparison.csv"),  row.names=FALSE)
  write.csv(params_slb, file.path(out_dir, "M1_Slb_htn16lb_all_params.csv"),  row.names=FALSE)
  saveRDS(msm_slb,      file.path(out_dir, "M1_msm_Slb.rds"))
  cat("Saved: M1_Slb_htn16lb_comparison.csv, M1_Slb_htn16lb_all_params.csv\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# S_ub: htn_baseline beta(1->6) fixed at upper 95% CI bound
# ═══════════════════════════════════════════════════════════════════════════════
cat(sprintf("\n\n====== S_ub: htn_baseline beta(1->6) fixed at upper CI = %.4f (HR=%.4f) ======\n",
            beta_ub, exp(beta_ub)))

# covinits: same as full model for all betas, except htn_baseline(1->6) = beta_ub
htn_covinits_ub <- htn_betas
htn_covinits_ub[htn16_pos] <- beta_ub

msm_sub <- tryCatch(
  msm(state ~ age, subject = patient_id, data = dat,
      qmatrix   = V_full_ref, covariates = COV_FULL, gen.inits = FALSE,
      covinits  = list(sex = sex_betas, htn_baseline = htn_covinits_ub),
      fixedpars = htn16_idx,
      control   = list(fnscale = 4000, maxit = 3000, reltol = 1e-8)),
  error = function(e) { cat("ERROR S_ub:", conditionMessage(e), "\n"); NULL }
)

params_sub <- NULL
if (!is.null(msm_sub)) {
  cat("S_ub converged.\n")
  params_sub <- extract_all(msm_sub, ref_cov, TRANS_Q, TRANS_NM, COV_NAMES)
  print_params(params_sub, sprintf("S_ub: htn_baseline beta(1->6) fixed = %.4f (upper CI)", beta_ub))
  cmp_sub <- compare_cov(params_full, params_sub, "Sub")
  cat("\nParameter changes (full -> S_ub):\n"); print_comparison(cmp_sub, "Sub")
  write.csv(cmp_sub,    file.path(out_dir, "M1_Sub_htn16ub_comparison.csv"),  row.names=FALSE)
  write.csv(params_sub, file.path(out_dir, "M1_Sub_htn16ub_all_params.csv"),  row.names=FALSE)
  saveRDS(msm_sub,      file.path(out_dir, "M1_msm_Sub.rds"))
  cat("Saved: M1_Sub_htn16ub_comparison.csv, M1_Sub_htn16ub_all_params.csv\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Side-by-side summary: Full / S_lb / S_ub
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n\n====== Side-by-Side: Full / S_lb / S_ub ======\n")
cat("  [*]=significant  [ ]=not significant  [F]=fixed\n\n")

get_row <- function(params, trans, cov) {
  if (is.null(params)) return(list(beta=NA, sig=NA, fix=NA))
  r <- params[params$block=="covariate" & params$transition==trans & params$covariate==cov, ]
  if (nrow(r)==0) return(list(beta=NA, sig=NA, fix=NA))
  list(beta=r$beta, sig=r$significant, fix=r$fixed)
}
fmt_cell <- function(x) {
  if (is.na(x$beta)) return(sprintf("%-20s", "—"))
  tag <- if (isTRUE(x$fix)) "[F]" else if (isTRUE(x$sig)) "[*]" else "[ ]"
  sprintf("%7.4f %-4s        ", x$beta, tag)
}

cv_full <- params_full[params_full$block=="covariate", ]
cat(sprintf("  %-8s %-14s  %-20s  %-20s  %-20s\n",
            "trans","covariate","Full",
            sprintf("S_lb(β=%.4f)", beta_lb),
            sprintf("S_ub(β=%.4f)", beta_ub)))
cat(sprintf("  %s\n", strrep("-", 84)))

sig_rows <- list()
for (i in seq_len(nrow(cv_full))) {
  tr <- cv_full$transition[i]; cv <- cv_full$covariate[i]
  r0  <- list(beta=cv_full$beta[i], sig=cv_full$significant[i], fix=cv_full$fixed[i])
  rlb <- get_row(params_slb, tr, cv)
  rub <- get_row(params_sub, tr, cv)
  cat(sprintf("  %-8s %-14s  %s  %s  %s\n",
              tr, cv, fmt_cell(r0), fmt_cell(rlb), fmt_cell(rub)))
  sig_rows[[i]] <- data.frame(
    transition=tr, covariate=cv,
    beta_full=round(r0$beta,4),  sig_full=r0$sig,  fixed_full=r0$fix,
    beta_Slb=round(rlb$beta,4),  sig_Slb=rlb$sig,  fixed_Slb=rlb$fix,
    beta_Sub=round(rub$beta,4),  sig_Sub=rub$sig,  fixed_Sub=rub$fix,
    stringsAsFactors=FALSE)
}
sig_df <- do.call(rbind, sig_rows)
write.csv(sig_df, file.path(out_dir, "M1_sensitivity_significance_table.csv"), row.names=FALSE)
cat(sprintf("\nSaved: M1_sensitivity_significance_table.csv\n"))
cat("\n===== Done =====\n")
sink()
