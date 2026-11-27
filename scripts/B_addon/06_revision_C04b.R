#!/usr/bin/env Rscript
# C04b: Outcome-conditional (Lash) probabilistic bias analysis for NOAF
# vs no-AF misclassification. Replaces C04 / B05 as a *corrected* estimate
# (those used outcome-independent relabeling, which pulls the HR toward 1
# and must not be cited as a correction). Does NOT overwrite C04 or B05 files.
#
# Writes ONLY to output_addon/ (and INTERNAL_ONLY therein).
# Never writes to output/ or results_final/. Never prints individual rows.
# READ_ONLY attach of catalog + mimic.duckdb. No new packages. Seed 20261120.
#
# Usage (project root):
#   nohup Rscript --vanilla scripts/B_addon/06_revision_C04b.R \
#        > logs/B_addon_revision_C04b.out 2>&1 &

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(survival)
})

# Paths come from environment variables / auto-detection (see scripts/00_helpers.R and README.md).
source(if (file.exists("scripts/00_helpers.R")) "scripts/00_helpers.R" else {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (!length(f)) stop("Set NOAF_REPO and run from the repository root (README.md).")
  file.path(dirname(dirname(dirname(normalizePath(f[1])))), "scripts/00_helpers.R")
}, local = FALSE)
addon_dir <- file.path(proj, "output_addon")
out_dir <- addon_dir
int_dir <- file.path(addon_dir, "INTERNAL_ONLY")
fig_dir <- file.path(addon_dir, "figures")
log_dir <- file.path(proj, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(log_dir, "B_addon_revision_C04b.log")
if (file.exists(logf)) file.remove(logf)
SEED <- 20261120L
set.seed(SEED)
log("START revision C04b (Lash outcome-conditional PBA); seed", SEED)

n_small <- function(x) { x <- as.numeric(x); !is.na(x) & x >= 1 & x <= 9 }
write_internal <- function(df, name) {
  write.csv(df, file.path(int_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
# Count mask: 1-9 -> <=9, 0 stays 0. Never mask design parameters (Se/Sp, seeds, B, horizons).
# Never mask infeasible_prop or PPV/NPV (these are not person-counts).
write_public <- function(df, name, count_cols = NULL) {
  pub <- df
  never_count <- c("horizon_days", "cut_lo", "cut_hi", "quantile_type",
                   "se", "sp", "Se", "Sp", "Se_death", "Se_surv", "Sp_death", "Sp_surv",
                   "seed", "seed_base", "seed_cell", "B", "B_ok", "scenario_id",
                   "youden_death", "youden_surv")
  if (is.null(count_cols)) {
    count_cols <- names(pub)[grepl("^(n_|n$)|_n$", names(pub)) | names(pub) %in% c("n", "N", "Freq", "n_event")]
  }
  count_cols <- setdiff(count_cols, never_count)
  small_by_col <- list()
  row_denom_small <- rep(FALSE, nrow(df))
  for (cc in count_cols) {
    if (!cc %in% names(df)) next
    sm <- n_small(df[[cc]])
    small_by_col[[cc]] <- sm
    if (cc %in% c("n", "N", "n_total", "n_event", "n_atrisk")) row_denom_small <- row_denom_small | sm
    x <- as.numeric(df[[cc]])
    pub[[cc]] <- ifelse(is.na(x), "",
                        ifelse(sm, "\u22649", format(as.integer(round(x)), scientific = FALSE, trim = TRUE)))
  }
  rcols <- names(df)[grepl("pct_|rate_|per_100|prop_|proportion|incidence", names(df), ignore.case = TRUE)]
  rcols <- setdiff(rcols, c("infeasible_prop", never_count))
  for (cc in rcols) {
    mask <- row_denom_small
    for (cand in names(small_by_col)) mask <- mask | small_by_col[[cand]]
    pub[[cc]] <- ifelse(mask, "\u22649", as.character(df[[cc]]))
  }
  write.csv(pub, file.path(out_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_both <- function(df, name, count_cols = NULL) {
  write_internal(df, name); write_public(df, name, count_cols)
}

# ---------------- Replay locked data build (no Flow / no official writes) ----------------
log("replay 01_run_main.R through stay-level analytic frame (in3p)")
main_src <- file.path(repo, "scripts/01_run_main.R")
src <- readLines(main_src, encoding = "UTF-8")
start_i <- grep("^catalog <- file.path", src)[1]
stop_i <- grep("# ---------------- Flow ----------------", src)[1]
if (is.na(start_i) || is.na(stop_i)) stop("could not locate catalog/Flow markers")
chunk <- src[start_i:(stop_i - 1)]
chunk <- chunk[!grepl("on\\.exit\\(.*dbDisconnect", chunk)]
eval(parse(text = paste(chunk, collapse = "\n")))
# Re-assert addon output paths in case the replay touched out_dir
out_dir <- addon_dir
int_dir <- file.path(addon_dir, "INTERNAL_ONLY")
log("analytic frame in memory; in3p defined; n=", sum(in3p), "DM", sum(in3p & dat$has_dm == 1))
try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
log("disconnected DuckDB (READ_ONLY session ended); C04b needs no further SQL")

# Counting-process constructor (same as locked / C04)
make_tvc <- function(df, time_var, status_var) {
  d <- df
  d$tstop0 <- d[[time_var]]
  d$ev <- d[[status_var]]
  d <- d[!is.na(d$tstop0) & d$tstop0 > 0, ]
  d$g <- as.character(d$af_group)
  split_noaf <- d$g == "NOAF" & !is.na(d$t_af_days) & d$t_af_days < d$tstop0
  a <- d
  a$tstart <- 0
  a$tstop <- ifelse(split_noaf, d$t_af_days, d$tstop0)
  a$event <- ifelse(split_noaf, 0L, d$ev)
  a$exp <- ifelse(d$g == "NOAF", "no_AF", d$g)
  b <- d[split_noaf, , drop = FALSE]
  if (nrow(b)) {
    b$tstart <- b$t_af_days
    b$tstop <- b$tstop0
    b$event <- b$ev
    b$exp <- "NOAF"
    out <- rbind(a, b)
  } else out <- a
  out$exp <- factor(out$exp, levels = c("no_AF", "NOAF", "preexisting"))
  out[out$tstop > out$tstart, ]
}
drop_noaf_missing_onset <- function(d, lab = "") {
  noaf_miss <- !is.na(d$af_group) & as.character(d$af_group) == "NOAF" & (is.na(d$t_af_days) | d$t_af_days <= 0)
  if (any(noaf_miss)) {
    log("drop NOAF missing t_af_days:", sum(noaf_miss), lab)
    d <- d[!noaf_miss, ]
  }
  d
}
restore_factors <- function(d) {
  if ("sex" %in% names(d) && !is.factor(d$sex)) d$sex <- factor(d$sex, levels = c("F", "M"))
  if ("race_coll" %in% names(d) && !is.factor(d$race_coll))
    d$race_coll <- factor(d$race_coll, levels = c("White", "Black", "Hispanic", "Asian", "Other or Unknown"))
  if ("adm_coll" %in% names(d) && !is.factor(d$adm_coll))
    d$adm_coll <- factor(d$adm_coll, levels = c("Emergency", "Elective", "Other"))
  if ("year_f" %in% names(d) && !is.factor(d$year_f)) d$year_f <- factor(d$year_f)
  if ("exp" %in% names(d) && !is.factor(d$exp)) d$exp <- factor(d$exp, levels = c("no_AF", "NOAF", "preexisting"))
  d
}
covars <- c("anchor_age", "sex", "race_coll", "adm_coll", "year_f", "sofa", "sapsii",
            "charlson_ex_dm", "sepsis3", "has_hf", "has_htn", "has_ckd", "cardiac_surgery",
            "vent24", "vaso24", "potassium_min", "mg_min", "gluc_mean")
form_rhs <- paste(c("exp", covars), collapse = " + ")

# ---------------- Lash / Fox outcome-conditional classification probabilities ----------------
# Within an outcome stratum, invert Se/Sp to expected true counts, then PPV/NPV.
# A = [n_obs_exp - (1-Sp)*n] / (Se+Sp-1);  PPV = Se*A / n_obs_exp;  NPV = Sp*B / n_obs_unexp.
# Infeasible if Youden <= 0, expected counts negative or > n, or PPV/NPV outside [0,1].
lash_stratum <- function(n_obs_exp, n_obs_unexp, Se, Sp, tol = 1e-10) {
  n <- n_obs_exp + n_obs_unexp
  out <- list(feasible = FALSE, reason = "ok", n = n,
              n_obs_exp = n_obs_exp, n_obs_unexp = n_obs_unexp,
              Se = Se, Sp = Sp, youden = Se + Sp - 1,
              prev_obs = NA_real_, E_true_exp = NA_real_, E_true_unexp = NA_real_,
              prev_true = NA_real_, PPV = NA_real_, NPV = NA_real_)
  if (is.na(Se) || is.na(Sp) || Se < 0 || Se > 1 || Sp < 0 || Sp > 1) {
    out$reason <- "se_sp_out"; return(out)
  }
  if (n <= 0) { out$reason <- "n_zero"; return(out) }
  if (out$youden <= tol) { out$reason <- "youden_le_0"; return(out) }
  out$prev_obs <- n_obs_exp / n
  A <- (n_obs_exp - (1 - Sp) * n) / out$youden
  B <- n - A
  if (abs(A) <= tol) A <- 0
  if (abs(B) <= tol) B <- 0
  if (A < -tol || B < -tol) {
    out$reason <- "neg_true_count"; out$E_true_exp <- A; out$E_true_unexp <- B
    return(out)
  }
  if (A > n + tol || B > n + tol) {
    out$reason <- "true_count_gt_n"; out$E_true_exp <- A; out$E_true_unexp <- B
    return(out)
  }
  A <- min(max(A, 0), n)
  B <- n - A
  out$E_true_exp <- A
  out$E_true_unexp <- B
  out$prev_true <- A / n
  PPV <- if (n_obs_exp > 0) (Se * A) / n_obs_exp else NA_real_
  NPV <- if (n_obs_unexp > 0) (Sp * B) / n_obs_unexp else NA_real_
  okp <- function(p) is.na(p) || (p >= -tol && p <= 1 + tol)
  if (!okp(PPV)) { out$reason <- "ppv_out"; out$PPV <- PPV; out$NPV <- NPV; return(out) }
  if (!okp(NPV)) { out$reason <- "npv_out"; out$PPV <- PPV; out$NPV <- NPV; return(out) }
  if (!is.na(PPV)) PPV <- min(max(PPV, 0), 1)
  if (!is.na(NPV)) NPV <- min(max(NPV, 0), 1)
  out$PPV <- PPV
  out$NPV <- NPV
  out$feasible <- TRUE
  out
}

# Formula self-check (aggregates only; no data)
.chk1 <- lash_stratum(200, 800, 0.80, 0.95)
.chk2 <- lash_stratum(10, 990, 0.80, 0.95)
if (!isTRUE(.chk1$feasible && abs(.chk1$E_true_exp - 200) < 1e-8 && abs(.chk1$PPV - 0.80) < 1e-8))
  stop("lash_stratum self-check 1 failed")
if (!isTRUE((!.chk2$feasible) && .chk2$reason == "neg_true_count"))
  stop("lash_stratum self-check 2 failed")
log("lash_stratum self-check ok")

# ---------------- Complete-case diabetes 28-day sample (same as C04) ----------------
cc_vars <- c("time_28", "status_28", "af_group", covars)
dcc <- dat[in3p & dat$has_dm == 1, ]
dcc <- dcc[stats::complete.cases(dcc[, cc_vars]), ]
dcc <- drop_noaf_missing_onset(dcc, "C04b")
dcc$af_group <- factor(as.character(dcc$af_group), levels = c("no_AF", "NOAF", "preexisting"))
dcc$status_28 <- as.integer(dcc$status_28)

n_cc <- nrow(dcc)
n_noaf_cc <- sum(dcc$af_group == "NOAF")
n_noaf_grp <- sum(dcc$af_group %in% c("no_AF", "NOAF"))
n_pre <- sum(dcc$af_group == "preexisting")
# Outcome-stratified 2x2 among reclassifiable stays (pre-existing AF excluded)
is_death <- dcc$status_28 == 1
is_noaf <- dcc$af_group == "NOAF"
is_no   <- dcc$af_group == "no_AF"
n_obs_NOAF_death <- sum(is_noaf & is_death)
n_obs_noAF_death <- sum(is_no & is_death)
n_obs_NOAF_surv  <- sum(is_noaf & !is_death)
n_obs_noAF_surv  <- sum(is_no & !is_death)
n_death_reclass <- n_obs_NOAF_death + n_obs_noAF_death
n_surv_reclass  <- n_obs_NOAF_surv + n_obs_noAF_surv
n_event_cc <- sum(is_death)

# Onset pools for false negatives, stratified by 28-day outcome (assumption; 未确认)
pool_death <- dcc$t_af_days[is_noaf & is_death & !is.na(dcc$t_af_days)]
pool_surv  <- dcc$t_af_days[is_noaf & !is_death & !is.na(dcc$t_af_days)]
n_pool_death <- length(pool_death)
n_pool_surv  <- length(pool_surv)
log("C04b CC n=", n_cc, "NOAF", n_noaf_cc, "preexisting", n_pre,
    "events", n_event_cc,
    "2x2 death NOAF/noAF", n_obs_NOAF_death, n_obs_noAF_death,
    "2x2 surv NOAF/noAF", n_obs_NOAF_surv, n_obs_noAF_surv,
    "onset pools death/surv", n_pool_death, n_pool_surv)

id_noaf_d <- which(is_noaf & is_death)
id_noaf_s <- which(is_noaf & !is_death)
id_no_d   <- which(is_no & is_death)
id_no_s   <- which(is_no & !is_death)

# Cox on one reclassified complete-case frame; return logHR and SE (no row printing)
fit_one_cc <- function(d) {
  long <- restore_factors(make_tvc(d, "time_28", "status_28"))
  long <- long[stats::complete.cases(long[, c("tstart", "tstop", "event", "exp", covars)]), ]
  if (nrow(long) < 50 || sum(long$event) < 10)
    return(c(loghr = NA_real_, se = NA_real_))
  f <- as.formula(paste("Surv(tstart, tstop, event) ~", form_rhs))
  fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(loghr = NA_real_, se = NA_real_))
  cf <- coef(fit)
  if (!"expNOAF" %in% names(cf)) return(c(loghr = NA_real_, se = NA_real_))
  vc <- try(vcov(fit), silent = TRUE)
  se <- if (!inherits(vc, "try-error") && "expNOAF" %in% colnames(vc))
    sqrt(vc["expNOAF", "expNOAF"]) else NA_real_
  c(loghr = unname(cf["expNOAF"]), se = unname(se))
}

obs_fit <- fit_one_cc(dcc)
hr_obs <- if (is.finite(obs_fit["loghr"])) exp(unname(obs_fit["loghr"])) else NA_real_
se_obs <- unname(obs_fit["se"])
log("C04b observed CC HR NOAF vs no_AF", hr_obs, "se_loghr", se_obs)

# Sample-level aggregates (written once; not overwritten onto C04 files)
samp <- data.frame(
  n = n_cc, n_NOAF = n_noaf_cc, n_no_AF = sum(is_no), n_preexisting = n_pre,
  n_noAF_plus_NOAF = n_noaf_grp, n_event = n_event_cc,
  n_obs_NOAF_death = n_obs_NOAF_death, n_obs_noAF_death = n_obs_noAF_death,
  n_obs_NOAF_surv = n_obs_NOAF_surv, n_obs_noAF_surv = n_obs_noAF_surv,
  n_death_reclass = n_death_reclass, n_surv_reclass = n_surv_reclass,
  n_onset_pool_death = n_pool_death, n_onset_pool_surv = n_pool_surv,
  hr_observed_cc = hr_obs, se_loghr_observed_cc = se_obs,
  stringsAsFactors = FALSE
)
write_both(samp, "C04b_sample_counts.csv",
           c("n", "n_NOAF", "n_no_AF", "n_preexisting", "n_noAF_plus_NOAF", "n_event",
             "n_obs_NOAF_death", "n_obs_noAF_death", "n_obs_NOAF_surv", "n_obs_noAF_surv",
             "n_death_reclass", "n_surv_reclass", "n_onset_pool_death", "n_onset_pool_surv"))

# ---------------- Scenario list ----------------
# Main grid: non-differential Se/Sp (same in death and survival).
# Differential Se (Sp remains non-differential, same Sp grid):
#   Se_death = Se_surv + 0.05, Se_surv in 0.60-0.90
#   Se_death = Se_surv + 0.10, Se_surv in 0.60-0.85
#   reverse: Se_death = Se_surv - 0.05, Se_surv in 0.65-0.95
se_grid <- seq(0.60, 0.95, by = 0.05)
sp_grid <- seq(0.95, 0.99, by = 0.01)
B_mc <- 200L

scen_spec <- list(
  list(id = 1L, name = "nondifferential",
       se_surv_vals = se_grid, d_se = 0,
       label = "Non-differential: Se_death = Se_surv, Sp_death = Sp_surv. Se 0.60-0.95 by 0.05, Sp 0.95-0.99 by 0.01."),
  list(id = 2L, name = "Se_death_plus_0.05",
       se_surv_vals = seq(0.60, 0.90, by = 0.05), d_se = 0.05,
       label = "Differential: Se_death = Se_surv + 0.05 (deaths ascertained more completely); Sp non-differential, same Sp grid. Se_surv 0.60-0.90."),
  list(id = 3L, name = "Se_death_plus_0.10",
       se_surv_vals = seq(0.60, 0.85, by = 0.05), d_se = 0.10,
       label = "Differential: Se_death = Se_surv + 0.10; Sp non-differential, same Sp grid. Se_surv 0.60-0.85."),
  list(id = 4L, name = "Se_death_minus_0.05",
       se_surv_vals = seq(0.65, 0.95, by = 0.05), d_se = -0.05,
       label = "Reverse differential: Se_death = Se_surv - 0.05 (deaths ascertained less completely); Sp non-differential, same Sp grid. Se_surv 0.65-0.95.")
)
scen_df <- do.call(rbind, lapply(scen_spec, function(s) {
  data.frame(scenario_id = s$id, scenario = s$name, delta_Se_death_minus_surv = s$d_se,
             Se_surv_min = min(s$se_surv_vals), Se_surv_max = max(s$se_surv_vals),
             Sp_min = min(sp_grid), Sp_max = max(sp_grid),
             n_Se_surv = length(s$se_surv_vals), n_Sp = length(sp_grid),
             n_cells = length(s$se_surv_vals) * length(sp_grid),
             description = s$label, stringsAsFactors = FALSE)
}))
# n_Se_surv / n_Sp / n_cells are grid sizes, not person counts — do not mask
write_both(scen_df, "C04b_scenario_list.csv", count_cols = character(0))
log("C04b scenarios:", nrow(scen_df), "total cells", sum(scen_df$n_cells), "B", B_mc)

# Reclassify one draw given stratum-specific PPV/NPV; onset for FN is outcome-stratified
reclassify_one <- function(PPV_d, NPV_d, PPV_s, NPV_s) {
  x <- dcc
  do_fp <- function(x, ids, PPV) {
    if (!length(ids) || is.na(PPV)) return(x)
    fp <- stats::runif(length(ids)) > PPV
    if (any(fp)) {
      x$af_group[ids[fp]] <- "no_AF"
      x$t_af_days[ids[fp]] <- NA_real_
    }
    x
  }
  do_fn <- function(x, ids, NPV, pool) {
    if (!length(ids) || is.na(NPV)) return(x)
    use_pool <- pool
    if (!length(use_pool)) use_pool <- c(pool_death, pool_surv)
    if (!length(use_pool)) return(x)
    fn <- stats::runif(length(ids)) > NPV
    if (any(fn)) {
      x$af_group[ids[fn]] <- "NOAF"
      ta <- sample(use_pool, sum(fn), replace = TRUE)
      tmax <- pmax(x$time_28[ids[fn]] - 1 / 48, 1 / 48)
      x$t_af_days[ids[fn]] <- pmin(ta, tmax)
    }
    x
  }
  x <- do_fp(x, id_noaf_d, PPV_d)
  x <- do_fp(x, id_noaf_s, PPV_s)
  x <- do_fn(x, id_no_d, NPV_d, pool_death)
  x <- do_fn(x, id_no_s, NPV_s, pool_surv)
  x$af_group <- factor(as.character(x$af_group), levels = c("no_AF", "NOAF", "preexisting"))
  x
}

run_cell <- function(scenario, scenario_id, Se_surv, Se_death, Sp_surv, Sp_death,
                     i_se, j_sp) {
  subseed <- as.integer(SEED + 100000L * scenario_id + 1000L * i_se + 10L * j_sp)
  ld <- lash_stratum(n_obs_NOAF_death, n_obs_noAF_death, Se_death, Sp_death)
  ls <- lash_stratum(n_obs_NOAF_surv,  n_obs_noAF_surv,  Se_surv,  Sp_surv)
  reason <- if (!ld$feasible && !ls$feasible) paste(ld$reason, ls$reason, sep = "|")
            else if (!ld$feasible) paste0("death:", ld$reason)
            else if (!ls$feasible) paste0("surv:", ls$reason)
            else "ok"
  feasible <- isTRUE(ld$feasible && ls$feasible)
  rec <- data.frame(
    scenario = scenario, scenario_id = scenario_id,
    Se_surv = Se_surv, Se_death = Se_death, Sp_surv = Sp_surv, Sp_death = Sp_death,
    cell_feasible = as.integer(feasible), infeasible_reason = reason,
    youden_death = ld$youden, youden_surv = ls$youden,
    prev_obs_death = ld$prev_obs, prev_obs_surv = ls$prev_obs,
    E_true_NOAF_death = ld$E_true_exp, E_true_NOAF_surv = ls$E_true_exp,
    prev_true_death = ld$prev_true, prev_true_surv = ls$prev_true,
    PPV_death = ld$PPV, NPV_death = ld$NPV,
    PPV_surv = ls$PPV, NPV_surv = ls$NPV,
    n = n_cc, n_NOAF = n_noaf_cc, n_noAF_plus_NOAF = n_noaf_grp, n_event = n_event_cc,
    n_obs_NOAF_death = n_obs_NOAF_death, n_obs_noAF_death = n_obs_noAF_death,
    n_obs_NOAF_surv = n_obs_NOAF_surv, n_obs_noAF_surv = n_obs_noAF_surv,
    B = B_mc, B_ok = 0L, n_infeasible = if (feasible) 0L else B_mc,
    n_cox_fail = 0L, infeasible_prop = if (feasible) 0 else 1,
    seed_base = SEED, seed_cell = subseed,
    hr_observed_cc = hr_obs,
    hr_sys_median = NA_real_, hr_sys_mean = NA_real_,
    hr_sys_lcl = NA_real_, hr_sys_ucl = NA_real_,
    hr_total_median = NA_real_, hr_total_mean = NA_real_,
    hr_total_lcl = NA_real_, hr_total_ucl = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!feasible) {
    log("C04b", scenario, "Se_s", Se_surv, "Se_d", Se_death, "Sp", Sp_death,
        "INFEASIBLE", reason)
    return(rec)
  }
  set.seed(subseed)
  hs <- rep(NA_real_, B_mc)
  ht <- rep(NA_real_, B_mc)
  n_cox_fail <- 0L
  for (b in seq_len(B_mc)) {
    x <- reclassify_one(ld$PPV, ld$NPV, ls$PPV, ls$NPV)
    ft <- fit_one_cc(x)
    if (!is.finite(ft["loghr"])) {
      n_cox_fail <- n_cox_fail + 1L
      next
    }
    hs[b] <- unname(ft["loghr"])
    if (is.finite(ft["se"]) && ft["se"] >= 0)
      ht[b] <- unname(ft["loghr"]) + stats::rnorm(1L, mean = 0, sd = ft["se"])
  }
  ok <- is.finite(hs)
  okt <- is.finite(ht)
  rec$B_ok <- sum(ok)
  rec$n_cox_fail <- n_cox_fail
  rec$n_infeasible <- 0L
  rec$infeasible_prop <- 0
  if (any(ok)) {
    rec$hr_sys_median <- exp(stats::median(hs[ok]))
    rec$hr_sys_mean <- exp(mean(hs[ok]))
    rec$hr_sys_lcl <- exp(unname(stats::quantile(hs[ok], 0.025, names = FALSE)))
    rec$hr_sys_ucl <- exp(unname(stats::quantile(hs[ok], 0.975, names = FALSE)))
  }
  if (any(okt)) {
    rec$hr_total_median <- exp(stats::median(ht[okt]))
    rec$hr_total_mean <- exp(mean(ht[okt]))
    rec$hr_total_lcl <- exp(unname(stats::quantile(ht[okt], 0.025, names = FALSE)))
    rec$hr_total_ucl <- exp(unname(stats::quantile(ht[okt], 0.975, names = FALSE)))
  }
  log("C04b", scenario, "Se_s", Se_surv, "Se_d", Se_death, "Sp", Sp_death,
      "sys median", if (any(ok)) sprintf("%.2f", rec$hr_sys_median) else "NA",
      "total", if (any(okt)) sprintf("%.2f (%.2f-%.2f)", rec$hr_total_median,
                                     rec$hr_total_lcl, rec$hr_total_ucl) else "NA",
      "B_ok", rec$B_ok)
  rec
}

# Map Se_surv to a stable 1-based index on the parent se_grid (0.60, 0.65, ..., 0.95)
se_index <- function(Se_surv) match(sprintf("%.2f", Se_surv), sprintf("%.2f", se_grid))
sp_index <- function(Sp) match(sprintf("%.2f", Sp), sprintf("%.2f", sp_grid))

set.seed(SEED)
rows <- list()
for (s in scen_spec) {
  for (Se_surv in s$se_surv_vals) {
    Se_death <- Se_surv + s$d_se
    i_se <- se_index(Se_surv)
    if (is.na(i_se)) stop("Se_surv not on parent grid: ", Se_surv)
    for (Sp in sp_grid) {
      j_sp <- sp_index(Sp)
      rows[[length(rows) + 1]] <- run_cell(
        scenario = s$name, scenario_id = s$id,
        Se_surv = Se_surv, Se_death = Se_death,
        Sp_surv = Sp, Sp_death = Sp,
        i_se = i_se, j_sp = j_sp
      )
    }
  }
}
c04b <- do.call(rbind, rows)
count_cols <- c("n", "n_NOAF", "n_noAF_plus_NOAF", "n_event",
                "n_obs_NOAF_death", "n_obs_noAF_death",
                "n_obs_NOAF_surv", "n_obs_noAF_surv",
                "n_infeasible", "n_cox_fail")
# E_true_* can be fractional; leave unrounded in the file, but if we pass them as
# count columns the public writer would integer-round. Keep them out of count_cols.
write_both(c04b, "C04b_misclassification_grid.csv", count_cols)
log("C04b grid written, cells", nrow(c04b),
    "feasible", sum(c04b$cell_feasible == 1),
    "infeasible", sum(c04b$cell_feasible == 0))

c04b_method <- data.frame(
  item = "C04b",
  method = paste(
    "Outcome-conditional quantitative bias analysis (Lash / Fox probabilistic bias analysis)",
    "for misclassification of NOAF vs no AF in the diabetes, 28-day, complete-case time-varying Cox model",
    "(same covariates as the locked main model). This supersedes C04 and B05 as a corrected-HR analysis:",
    "C04 and B05 reclassified records with a single PPV/NPV computed from overall prevalence, independent of the outcome,",
    "which is equivalent to random relabeling and pulls the HR toward 1; those files are retained but should not be cited as corrected estimates.",
    "Pre-existing AF is held fixed and is not reclassified, and is excluded from the outcome-stratified 2x2.",
    "Within each 28-day outcome stratum (death / survival), expected true NOAF counts are obtained by inverting Se and Sp:",
    "E_true = [n_obs_NOAF - (1-Sp)*n] / (Se+Sp-1). Then PPV = Se*E_true / n_obs_NOAF and NPV = Sp*(n-E_true) / n_obs_noAF.",
    "A cell is infeasible if either stratum yields a negative expected count, an expected count outside [0,n], Youden <= 0, or PPV/NPV outside [0,1];",
    "infeasible_prop is 1 for such cells (no Cox fits) and 0 otherwise (Se/Sp are fixed, so infeasibility is deterministic given the observed 2x2).",
    "n_cox_fail counts iterations in feasible cells where the Cox model did not return a finite NOAF log HR.",
    "Each feasible iteration: classified NOAF are recoded to no AF with probability 1-PPV_stratum (onset cleared);",
    "classified no AF are recoded to NOAF with probability 1-NPV_stratum.",
    "Onset times for false negatives are sampled with replacement from observed NOAF t_af_days in the *same 28-day outcome stratum*,",
    "capped strictly below follow-up end (time_28 - 1/48 day, or 1/48 if follow-up is shorter). If a stratum pool is empty, the combined pool is used.",
    "The sampled-onset rule is an assumption; true missed-onset times are 未确认.",
    "B=200 iterations per cell. Systematic-error interval: median and 2.5th-97.5th percentiles of exp(log HR) from the reclassified Cox fits.",
    "Total-uncertainty interval: for each successful fit, add a N(0, SE^2) draw to that iteration's log HR (SE from the reclassified Cox model) and take the same percentiles;",
    "this combines classification uncertainty with conventional sampling error (Lash / Fox). It is not a bootstrap of the cohort.",
    "Base seed 20261120; cell seed = 20261120 + 100000*scenario_id + 1000*i_Se + 10*j_Sp,",
    "where i_Se and j_Sp are 1-based indices of Se_surv and Sp on the parent grids (Se 0.60-0.95 by 0.05; Sp 0.95-0.99 by 0.01).",
    "scenario_id 1 = nondifferential; 2 = Se_death = Se_surv+0.05; 3 = Se_death = Se_surv+0.10; 4 = Se_death = Se_surv-0.05 (reverse).",
    "In differential scenarios Sp remains non-differential. This is not a multiple-imputation pooled HR."
  ),
  stringsAsFactors = FALSE
)
write_both(c04b_method, "C04b_method_note.csv")

sess <- data.frame(
  item = c("date", "timezone", "R", "survival", "duckdb", "seed", "B",
           "cohort", "method", "n_cells_planned"),
  value = c(as.character(Sys.time()), "Pacific/Apia (UTC+13)",
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(utils::packageVersion("survival")),
            as.character(utils::packageVersion("duckdb")),
            as.character(SEED), as.character(B_mc),
            "DM, 2008-2019, three-group in3p, complete-case; outcome 28-day death; contrast NOAF vs no AF, TVC Cox",
            "Lash outcome-conditional PBA; systematic + total (add N(0,SE^2)) intervals",
            as.character(sum(scen_df$n_cells))),
  stringsAsFactors = FALSE
)
write_both(sess, "C04b_session_info.csv")

of <- file.path(proj, "output", "01_flow_counts.csv")
if (file.exists(of)) log("official output/01_flow_counts.csv mtime still", as.character(file.info(of)$mtime))
# Confirm we did not touch old C04 public files
c04f <- file.path(out_dir, "C04_misclassification_grid.csv")
if (file.exists(c04f)) log("old C04 grid mtime unchanged check:", as.character(file.info(c04f)$mtime))
log("DONE revision C04b")
