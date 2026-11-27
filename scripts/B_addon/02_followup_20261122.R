#!/usr/bin/env Rscript
# CVD add-on follow-up 2026-11-22.
# Writes ONLY to output_addon/ (and INTERNAL_ONLY therein).
# Never writes to output/ or results_final/. Never prints individual rows.
# Does not install packages. Does not edit manuscripts.
#
# 1. Reconcile 48h landmark original-spec n/HR vs locked 07_hr_landmark48.csv
# 2. Re-export Fig 2 / Fig S1 Simon-Makuch (and Fig S2 CIF) without in-figure captions
# 3. Record Fig S2 CIF grouping (baseline vs time-varying)
# 4. Table S8 tertile cut-points
# 5. strong::ok_mult_cox with aligned num_covar; cox.zph on the same time-fixed model

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(survival); library(ggplot2)
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
logf <- file.path(log_dir, "B_addon_followup.log")
if (file.exists(logf)) file.remove(logf)
SEED <- 20261120L
set.seed(SEED)
log("START B_addon follow-up 2026-11-22; seed", SEED)

has_strong <- isTRUE(requireNamespace("strong", quietly = TRUE))
log("strong available:", has_strong, if (has_strong) as.character(packageVersion("strong")) else "")

n_small <- function(x) { x <- as.numeric(x); !is.na(x) & x <= 9 }
write_internal <- function(df, name) {
  write.csv(df, file.path(int_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_public <- function(df, name, count_cols = NULL) {
  pub <- df
  if (is.null(count_cols)) {
    count_cols <- names(pub)[grepl("^(n_|n$)|_n$", names(pub)) | names(pub) %in% c("n", "N", "Freq", "n_event")]
  }
  small_by_col <- list()
  row_denom_small <- rep(FALSE, nrow(df))
  for (cc in count_cols) {
    if (!cc %in% names(df)) next
    sm <- n_small(df[[cc]])
    small_by_col[[cc]] <- sm
    if (cc %in% c("n", "N", "n_total", "n_event", "n_atrisk", "n_remaining", "n_excluded")) {
      row_denom_small <- row_denom_small | sm
    }
    x <- as.numeric(df[[cc]])
    pub[[cc]] <- ifelse(is.na(x), "", ifelse(sm, "\u22649", format(as.integer(round(x)), scientific = FALSE, trim = TRUE)))
  }
  rcols <- names(df)[grepl("pct_|rate_|per_100|incidence", names(df), ignore.case = TRUE)]
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

# ---------------- Replay locked data build through stay-level frame (no official writes) ----------------
log("replay 01_run_main.R through stay-level analytic frame (in3p); no Flow/write sections")
main_src <- file.path(repo, "scripts/01_run_main.R")
src <- readLines(main_src, encoding = "UTF-8")
start_i <- grep("^catalog <- file.path", src)[1]
stop_i <- grep("# ---------------- Flow ----------------", src)[1]
if (is.na(start_i) || is.na(stop_i)) stop("could not locate catalog/Flow markers in 01_run_main.R")
chunk <- src[start_i:(stop_i - 1)]
chunk <- chunk[!grepl("on\\.exit\\(.*dbDisconnect", chunk)]
eval(parse(text = paste(chunk, collapse = "\n")))
log("analytic frame in memory; n first-stay rows not printed; in3p defined")
san_n <- sum(in3p)
san_dm <- sum(in3p & dat$has_dm == 1)
log("prep sanity three-group 2008-19 n=", san_n, "DM", san_dm)

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

covars <- c("anchor_age", "sex", "race_coll", "adm_coll", "year_f", "sofa", "sapsii",
            "charlson_ex_dm", "sepsis3", "has_hf", "has_htn", "has_ckd", "cardiac_surgery",
            "vent24", "vaso24", "potassium_min", "mg_min", "gluc_mean")
covars_num_intended <- c("anchor_age", "sofa", "sapsii", "charlson_ex_dm", "potassium_min", "mg_min", "gluc_mean")
bin_covars <- c("sepsis3", "has_hf", "has_htn", "has_ckd", "cardiac_surgery", "vent24", "vaso24")

drop_noaf_missing_onset <- function(d, lab = "") {
  noaf_miss <- !is.na(d$af_group) & as.character(d$af_group) == "NOAF" & (is.na(d$t_af_days) | d$t_af_days <= 0)
  if (any(noaf_miss)) {
    log("drop NOAF missing t_af_days:", sum(noaf_miss), lab)
    d <- d[!noaf_miss, ]
  }
  d
}

extract_hr <- function(fit, model, outcome, subset_lab, n = NULL, n_event = NULL) {
  sm <- summary(fit)$conf.int
  if (is.null(sm)) return(NULL)
  rn <- rownames(sm)
  keep <- grepl("^exp", rn)
  if (!any(keep)) return(NULL)
  data.frame(
    model = model, outcome = outcome, subset = subset_lab, term = sub("^exp", "", rn[keep]),
    hr = sm[keep, "exp(coef)"],
    lcl = sm[keep, "lower .95"],
    ucl = sm[keep, "upper .95"],
    p = summary(fit)$coefficients[keep, "Pr(>|z|)"],
    n = if (is.null(n)) fit$n else n,
    n_event = if (is.null(n_event)) as.integer(fit$nevent) else as.integer(n_event),
    stringsAsFactors = FALSE
  )
}

grp_n <- function(d) {
  g <- as.character(d$af_group)
  c(n_no_AF = sum(g == "no_AF"), n_NOAF = sum(g == "NOAF"), n_pre = sum(g == "preexisting"))
}
lm_grp_n <- function(d) {
  g <- as.character(d$lm_group)
  c(n_no_AF = sum(g == "no_AF"), n_NOAF = sum(g == "NOAF"), n_pre = sum(g == "preexisting"))
}

# =====================================================================
# Item 1: 48h landmark original-spec reconciliation (diabetes, 28-day)
# =====================================================================
log("Item 1: 48h landmark original-spec reconciliation")
d0 <- dat[in3p & dat$has_dm == 1, ]
n0 <- nrow(d0)
g0 <- grp_n(d0)
n_ev0 <- sum(d0$status_28 == 1, na.rm = TRUE)

# Locked / SAP: alive (followed) beyond 48h on the 28-day clock; no in-hospital requirement
d_alive <- d0[d0$time_28 > 2, ]
n_alive <- nrow(d_alive)
n_die48 <- n0 - n_alive
# Redundant second filter in locked code
d_alive2 <- d_alive[!(d_alive$status_28 == 1 & d_alive$time_28 <= 2), ]
n_redund <- nrow(d_alive) - nrow(d_alive2)

# Addon extra filter applied to 28-day/1-year old_style: still in hospital at 48h
# time_hosp = days from ICU intime date to hospital dischtime
still_hosp <- is.na(d_alive2$time_hosp) | d_alive2$time_hosp > 2
d_addon <- d_alive2[still_hosp, ]
d_disch <- d_alive2[!still_hosp, ]
n_extra <- nrow(d_alive2) - nrow(d_addon)
n_extra_ev <- sum(d_disch$status_28 == 1, na.rm = TRUE)

recode_late <- function(d) {
  d$lm_group <- as.character(d$af_group)
  d$lm_group[d$af_group == "NOAF" & (is.na(d$t_af_days) | d$t_af_days > 2)] <- "no_AF"
  d$lm_group <- factor(d$lm_group, levels = c("no_AF", "NOAF", "preexisting"))
  d$lm_time28 <- pmax(d$time_28 - 2, 1 / 24)
  d
}
d_lock <- recode_late(d_alive2)
d_add  <- recode_late(d_addon)
g_lock <- lm_grp_n(d_lock)
g_add  <- lm_grp_n(d_add)

cc_vars <- c("lm_time28", "status_28", "lm_group", covars)
d_lock_cc <- d_lock[stats::complete.cases(d_lock[, cc_vars]), ]
d_add_cc  <- d_add[stats::complete.cases(d_add[, cc_vars]), ]
n_lock_cc <- nrow(d_lock_cc)
n_add_cc  <- nrow(d_add_cc)
# Among the extra-excluded (alive at 48h, discharged by 48h), how many were complete-case in the locked sample?
d_disch_r <- recode_late(d_disch)
n_extra_in_lock_cc <- sum(stats::complete.cases(d_disch_r[, cc_vars]))

step_row <- function(step, spec, desc, n, n_excl, n_event, gn, note) {
  data.frame(step = step, spec = spec, description = desc,
             n = as.integer(n), n_excluded_this_step = as.integer(n_excl),
             n_event_28 = as.integer(n_event),
             n_no_AF = as.integer(gn["n_no_AF"]), n_NOAF = as.integer(gn["n_NOAF"]),
             n_pre = as.integer(gn["n_pre"]), note = note, stringsAsFactors = FALSE)
}

recon <- rbind(
  step_row(1, "common", "Diabetes, three-group, 2008-2019 (in3p)", n0, 0, n_ev0, g0,
           "Locked analytic frame; matches B00 and the locked cohort (n=14,329)."),
  step_row(2, "locked_and_SAP", "Alive at 48 h: time_28 > 2 days", n_alive, n_die48,
           sum(d_alive$status_28 == 1, na.rm = TRUE), grp_n(d_alive),
           "SAP v1: analysis set is those still alive 48 h after ICU intime. time_28 is 28-day all-cause follow-up (pmin of days_to_death or 28). time_28<=2 is essentially death within 48 h, not hospital discharge."),
  step_row(3, "locked_and_SAP", "Redundant filter: drop status_28==1 and time_28<=2",
           nrow(d_alive2), n_redund, sum(d_alive2$status_28 == 1, na.rm = TRUE), grp_n(d_alive2),
           "Already implied by time_28>2; should exclude 0 additional stays. Present in both locked 01_run_main.R and the add-on."),
  step_row(4, "addon_old_style_only", "EXTRA: still in hospital at 48 h (time_hosp > 2 or time_hosp missing)",
           nrow(d_addon), n_extra, sum(d_addon$status_28 == 1, na.rm = TRUE), grp_n(d_addon),
           paste0("This filter is NOT in SAP v1 or locked 01_run_main.R for 28-day/1-year landmark. time_hosp = date_diff(day, ICU intime, hospital dischtime). ",
                  "The add-on copied a 'still in hospital' rule from the 24 h / hospital-death landmark comment. ",
                  "n discharged by 48 h but still alive at 48 h = ", n_extra, " (of whom ", n_extra_ev, " later died by day 28).")),
  step_row(5, "locked_and_SAP", "Recode late NOAF (t_af>2 or missing onset) as no_AF; do not drop them",
           nrow(d_lock), 0, sum(d_lock$status_28 == 1, na.rm = TRUE), g_lock,
           "SAP v1 and locked code: NOAF only if t_af<=48 h; later AF remains no AF at the landmark. Sample size unchanged. Locked 07_landmark48_counts.csv is this step (before complete-case)."),
  step_row(5, "addon_old_style_only", "Same late-NOAF recode, after the extra hospital filter",
           nrow(d_add), 0, sum(d_add$status_28 == 1, na.rm = TRUE), g_add,
           "Same recode as locked; applied to the smaller in-hospital-at-48h sample."),
  step_row(6, "locked_and_SAP", "Complete-case on lm_time, status_28, lm_group, and the locked covars",
           n_lock_cc, nrow(d_lock) - n_lock_cc, sum(d_lock_cc$status_28 == 1, na.rm = TRUE), lm_grp_n(d_lock_cc),
           "Locked 07_hr_landmark48.csv: n=12,991, 1,354 deaths, NOAF HR 1.29 (1.02-1.63). Both locked and add-on are complete-case (not MI). Follow-up starts at day 2 (lm_time = pmax(time_28-2, 1/24))."),
  step_row(6, "addon_old_style_only", "Complete-case on the same covariates, after the extra hospital filter",
           n_add_cc, nrow(d_add) - n_add_cc, sum(d_add_cc$status_28 == 1, na.rm = TRUE), lm_grp_n(d_add_cc),
           paste0("Add-on B04 landmark_2d_old_style_cc reported n=11,781, HR 1.32. Difference vs locked complete-case n is ",
                  n_lock_cc - n_add_cc, "; of the ", n_extra, " stays dropped only by the extra hospital filter, ",
                  n_extra_in_lock_cc, " were complete-case in the locked model and account for this n gap."))
)
recon$verdict <- ""
recon$verdict[recon$spec == "locked_and_SAP" & recon$step == 6] <-
  "CORRECT per SAP v1 and locked 01_run_main.R. Use this specification (n should match 12,991)."
recon$verdict[recon$spec == "addon_old_style_only" & recon$step == 6] <-
  "INCORRECT for a rerun of the original locked 48 h landmark: the extra still-in-hospital filter is not in SAP or the locked code for 28-day/1-year outcomes."
write_both(recon, "B04_landmark48_reconciliation.csv",
           c("n", "n_excluded_this_step", "n_event_28", "n_no_AF", "n_NOAF", "n_pre"))

# Refit both specs to confirm HRs (complete-case, same formula as locked)
fit_lm_cc <- function(dcc, model) {
  dcc$exp <- dcc$lm_group
  f <- as.formula(paste("Surv(lm_time28, status_28) ~ exp +", paste(covars, collapse = " + ")))
  fit <- coxph(f, data = dcc, ties = "efron")
  extract_hr(fit, model, "death_28d", "DM", n = nrow(dcc), n_event = sum(dcc$status_28))
}
hr_lock <- fit_lm_cc(d_lock_cc, "landmark_48h_locked_spec_cc")
hr_add  <- fit_lm_cc(d_add_cc, "landmark_48h_addon_old_style_cc")
hr_both <- rbind(hr_lock, hr_add)
hr_both$locked_published_n <- 12991L
hr_both$locked_published_hr <- 1.28813618153825
hr_both$note <- ifelse(hr_both$model == "landmark_48h_locked_spec_cc",
                       "Should reproduce locked 07_hr_landmark48.csv (DM, 28-day).",
                       "Reproduces add-on B04 old_style; differs because of the extra hospital-stay filter.")
write_both(hr_both, "B04_landmark48_spec_rerun.csv", c("n", "n_event", "locked_published_n"))
log("Item 1 done; locked-spec CC n=", n_lock_cc, "addon-spec CC n=", n_add_cc)

# Keep locked-spec complete-case frame for item 5 (do not print)
lm48_dm_cc <- d_lock_cc

# =====================================================================
# Item 4: Table S8 tertile cut-points (before figures so a failure there cannot skip this)
# =====================================================================
log("Item 4: Table S8 glucose tertile cut-points")
dm28 <- dat[in3p & dat$has_dm == 1, ]
qs_cv <- as.numeric(stats::quantile(dm28$gluc_cv, c(1/3, 2/3), na.rm = TRUE, type = 7))
qs_mn <- as.numeric(stats::quantile(dm28$gluc_mean, c(1/3, 2/3), na.rm = TRUE, type = 7))
log("gluc_cv tertile cuts", qs_cv[1], qs_cv[2], "gluc_mean cuts", qs_mn[1], qs_mn[2])

ter_counts <- function(x, cuts, labels, varname, unit, cut_note) {
  xx <- x
  ter <- cut(xx, breaks = c(-Inf, cuts, Inf), labels = labels, right = TRUE, include.lowest = TRUE)
  # n with the variable observed (the population in which cuts were computed)
  tab_obs <- table(ter, useNA = "no")
  # complete-case n as in locked fit_tvc_cox (covars + time_28/status_28/af_group, then drop_noaf_missing_onset)
  d <- dm28
  d$ter <- ter
  rows <- list()
  for (lv in labels) {
    s <- d[!is.na(d$ter) & d$ter == lv, ]
    n_obs <- nrow(s)
    need <- c("time_28", "status_28", "af_group", covars)
    scc <- s[stats::complete.cases(s[, need]), ]
    scc <- drop_noaf_missing_onset(scc, paste("S8", varname, lv))
    rows[[length(rows) + 1]] <- data.frame(
      variable = varname,
      unit = unit,
      tertile = lv,
      cut_lo = if (lv == labels[1]) NA_real_ else cuts[1],
      cut_hi = if (lv == labels[3]) NA_real_ else if (lv == labels[1]) cuts[1] else cuts[2],
      interval = if (lv == labels[1]) sprintf("(-Inf, %s]", format(cuts[1], digits = 8, scientific = FALSE))
                 else if (lv == labels[2]) sprintf("(%s, %s]", format(cuts[1], digits = 8, scientific = FALSE),
                                                   format(cuts[2], digits = 8, scientific = FALSE))
                 else sprintf("(%s, Inf)", format(cuts[2], digits = 8, scientific = FALSE)),
      quantile_type = 7L,
      n_with_variable = as.integer(n_obs),
      n_complete_case = as.integer(nrow(scc)),
      n_event = as.integer(sum(scc$status_28 == 1, na.rm = TRUE)),
      n_missing_variable_in_DM_in3p = as.integer(sum(is.na(x))),
      cut_note = cut_note,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

s8 <- rbind(
  ter_counts(dm28$gluc_cv, qs_cv, c("CVt1_low", "CVt2", "CVt3_high"),
             "gluc_cv", "unitless (SD/mean of first-24h glucose)",
             "Cuts = quantile(..., c(1/3,2/3), na.rm=TRUE, type=7) among diabetes, 2008-19 three-group stays with non-missing gluc_cv, BEFORE complete-case. Locked 01_run_main.R lines ~1305-1310. Table S8 n is n_complete_case (after covars complete-case and drop of NOAF with missing onset)."),
  ter_counts(dm28$gluc_mean, qs_mn, c("meant1_low", "meant2", "meant3_high"),
             "gluc_mean", "mg/dL (first-24h mean glucose)",
             "Cuts = quantile(..., c(1/3,2/3), na.rm=TRUE, type=7) among diabetes, 2008-19 three-group stays with non-missing gluc_mean, BEFORE complete-case. Locked 01_run_main.R lines ~1306-1317. Table S8 n is n_complete_case.")
)
# Non-tertile rows for completeness (Table S8 also has hypo and SHR)
need <- c("time_28", "status_28", "af_group", covars)
hypo_rows <- lapply(c(0, 1), function(hv) {
  s <- dm28[dm28$hypo70 == hv, ]
  scc <- s[stats::complete.cases(s[, need]), ]
  scc <- drop_noaf_missing_onset(scc, paste("S8 hypo", hv))
  data.frame(variable = "hypo70", unit = "binary; any first-24h glucose <70 mg/dL",
             tertile = paste0("hypo70=", hv), cut_lo = NA_real_, cut_hi = NA_real_,
             interval = if (hv == 0) "no hypoglycaemia" else "any hypoglycaemia (<70 mg/dL)",
             quantile_type = NA_integer_,
             n_with_variable = as.integer(nrow(s)),
             n_complete_case = as.integer(nrow(scc)),
             n_event = as.integer(sum(scc$status_28 == 1, na.rm = TRUE)),
             n_missing_variable_in_DM_in3p = as.integer(sum(is.na(dm28$hypo70))),
             cut_note = "Not a tertile; included so Table S8 row n can be checked.",
             stringsAsFactors = FALSE)
})
shr_s <- dm28[!is.na(dm28$shr), ]
shr_cc <- shr_s[stats::complete.cases(shr_s[, need]), ]
shr_cc <- drop_noaf_missing_onset(shr_cc, "S8 SHR")
shr_row <- data.frame(variable = "shr", unit = "glucose / (28.7*HbA1c-46.7); subset with HbA1c",
                      tertile = "SHR subset", cut_lo = NA_real_, cut_hi = NA_real_,
                      interval = "has HbA1c and valid SHR (not tertiles)",
                      quantile_type = NA_integer_,
                      n_with_variable = as.integer(nrow(shr_s)),
                      n_complete_case = as.integer(nrow(shr_cc)),
                      n_event = as.integer(sum(shr_cc$status_28 == 1, na.rm = TRUE)),
                      n_missing_variable_in_DM_in3p = as.integer(sum(is.na(dm28$shr))),
                      cut_note = "Not a tertile; exploratory SHR subset as in locked 10_hr_dm_glucose_strata.csv.",
                      stringsAsFactors = FALSE)
s8 <- rbind(s8, do.call(rbind, hypo_rows), shr_row)
write_both(s8, "B11_tableS8_tertile_cutpoints.csv",
           c("n_with_variable", "n_complete_case", "n_event", "n_missing_variable_in_DM_in3p"))
log("Item 4 (S8 cuts) done")

# =====================================================================
# Item 3: CIF grouping note
# =====================================================================
cif_note <- data.frame(
  figure = "Fig. S2 (FigureS2_CIF_hospital_death)",
  grouping = "baseline_af_group",
  grouping_precise = "Baseline (time-fixed) af_group at ICU admission: no_AF / NOAF / pre-existing AF. Patients later classified as NOAF are placed in the NOAF curve from time zero, including the period before onset. The CIF is estimated with cmprsk::cuminc(ftime, fstatus, group) on hospital follow-up (status 1 = hospital death, 2 = live discharge as competing event, displayed to day 28), faceted by diabetes. It is NOT a Simon-Makuch / counting-process CIF and does not move patients from no AF to NOAF at onset.",
  implications = "Early hospital-death CIF in the NOAF panel can include time before AF onset (exposure is not time-varying). This is a descriptive baseline-group competing-risk plot; adjusted inference in B03 uses time-varying Fine-Gray and cause-specific Cox. Do not describe Fig. S2 as a time-varying exposure CIF. Contrast with Fig. 2 / Fig. S1, which do move patients at onset.",
  source_function = "plot_cif() in scripts/B_addon/01_run_addon.R (comments: 'Descriptive CIF by baseline AF group')",
  stringsAsFactors = FALSE
)
write_both(cif_note, "B03_cif_grouping_note.csv")
log("Item 3: CIF grouping is baseline af_group, not time-varying")

# =====================================================================
# Item 2 and S2: re-export figures without in-image captions
# =====================================================================
theme_clean <- function(base = 11) {
  theme_bw(base_size = base, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(8, 12, 8, 10),
      axis.title = element_text(size = base),
      axis.text = element_text(size = base - 1)
    )
}

save_fig_set <- function(p, name, width = 8.5, height = 5.0) {
  p <- p + labs(title = NULL, subtitle = NULL, caption = NULL)
  pdfp <- file.path(fig_dir, paste0(name, ".pdf"))
  pngp <- file.path(fig_dir, paste0(name, ".png"))
  tifp <- file.path(fig_dir, paste0(name, ".tiff"))
  ggsave(pdfp, p, width = width, height = height, device = grDevices::pdf, useDingbats = FALSE, bg = "white")
  ggsave(pngp, p, width = width, height = height, dpi = 300, device = "png", bg = "white")
  ggsave(tifp, p, width = width, height = height, dpi = 300, device = "tiff",
         compression = "lzw", bg = "white")
  log("wrote", name, "pdf/png/tiff")
}

pal <- c("No AF" = "#0072B2", "NOAF" = "#D55E00", "Pre-existing AF" = "#009E73")

sm_plot <- function(time_var, status_var, xmax, xlab, ylab, fname, fname_alias = NULL) {
  d <- dat[in3p, ]
  d$dmf <- factor(d$has_dm, levels = c(1, 0), labels = c("Diabetes", "No diabetes"))
  pieces <- list()
  for (dmv in levels(d$dmf)) {
    s <- d[d$dmf == dmv, ]
    s <- drop_noaf_missing_onset(s, paste("SM", dmv))
    long <- make_tvc(s, time_var, status_var)
    long$exp <- factor(long$exp, levels = c("no_AF", "NOAF", "preexisting"),
                       labels = c("No AF", "NOAF", "Pre-existing AF"))
    sf <- survfit(Surv(tstart, tstop, event) ~ exp, data = long)
    sm <- summary(sf)
    st <- if (is.null(sm$strata)) "all" else sub("^exp=", "", as.character(sm$strata))
    pieces[[length(pieces) + 1]] <- data.frame(dm = dmv, group = st, time = sm$time, surv = sm$surv,
                                               stringsAsFactors = FALSE)
  }
  pdat <- do.call(rbind, pieces)
  pdat$group <- factor(pdat$group, levels = c("No AF", "NOAF", "Pre-existing AF"))
  # no in-figure caption/notes; no numbers-at-risk were present previously (not added)
  p <- ggplot(pdat, aes(x = time, y = 1 - surv, color = group)) +
    geom_step(linewidth = 0.7) +
    facet_wrap(~dm) +
    scale_color_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = pretty(c(0, xmax), 6), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0, NA), clip = "on") +
    labs(x = xlab, y = ylab, title = NULL, subtitle = NULL, caption = NULL) +
    theme_clean()
  save_fig_set(p, fname)
  if (!is.null(fname_alias) && !identical(fname_alias, fname)) save_fig_set(p, fname_alias)
  invisible(p)
}

log("Item 2: re-export Simon-Makuch Fig 2 and Fig S1, no captions")
try(sm_plot("time_28", "status_28", 28,
            "Days from ICU admission", "Cumulative incidence of death",
            "Figure2_SimonMakuch_28day"), silent = FALSE)
try(sm_plot("time_1y", "status_1y", 365,
            "Days from ICU admission", "Cumulative incidence of death",
            "FigureS1_SimonMakuch_1year",
            fname_alias = "Figure2_SimonMakuch_1year"), silent = FALSE)

log("Also re-export Fig S2 CIF without in-figure caption (same baseline grouping)")
plot_cif <- function() {
  d <- dat[in3p, ]
  d$group <- factor(as.character(d$af_group), levels = c("no_AF", "NOAF", "preexisting"),
                    labels = c("No AF", "NOAF", "Pre-existing AF"))
  d$dmf <- factor(d$has_dm, levels = c(1, 0), labels = c("Diabetes", "No diabetes"))
  d$ftime <- pmax(d$time_hosp, 1 / 24)
  d$fstatus <- ifelse(d$status_hosp == 1, 1, 2)
  d$ftime28 <- pmin(d$ftime, 28)
  d$fstatus28 <- ifelse(d$ftime > 28, 0, d$fstatus)
  rows <- list()
  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    log("cmprsk not available; skip CIF redraw")
    return(NULL)
  }
  for (dmv in levels(d$dmf)) {
    s <- d[d$dmf == dmv, ]
    ci <- try(cmprsk::cuminc(s$ftime28, s$fstatus28, s$group, cencode = 0), silent = TRUE)
    if (inherits(ci, "try-error")) next
    nms <- names(ci)
    nms <- nms[nms != "Tests"]
    for (nm in nms) {
      if (!grepl(" 1$", nm)) next
      gp <- sub(" 1$", "", nm)
      xx <- ci[[nm]]
      rows[[length(rows) + 1]] <- data.frame(dm = dmv, group = gp, time = xx$time, cif = xx$est,
                                             stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  pdat <- do.call(rbind, rows)
  pdat$group <- factor(pdat$group, levels = c("No AF", "NOAF", "Pre-existing AF"))
  ggplot(pdat, aes(x = time, y = cif, color = group)) +
    geom_step(linewidth = 0.7) +
    facet_wrap(~dm) +
    scale_color_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = seq(0, 28, 7), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = function(z) sprintf("%.2f", z), expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(xlim = c(0, 28), ylim = c(0, NA), clip = "on") +
    labs(x = "Days from ICU admission", y = "Cumulative incidence of hospital death",
         title = NULL, subtitle = NULL, caption = NULL) +
    theme_clean()
}
p_cif <- try(plot_cif(), silent = TRUE)
if (!inherits(p_cif, "try-error") && !is.null(p_cif)) {
  save_fig_set(p_cif, "FigureS2_CIF_hospital_death")
} else log("CIF redraw failed", if (inherits(p_cif, "try-error")) as.character(p_cif)[1] else "null")

# =====================================================================
# Item 5: strong::ok_mult_cox on locked-spec 48h landmark (time-fixed), plus cox.zph
# =====================================================================
log("Item 5: strong::ok_mult_cox with aligned num_covar; cox.zph on same model")
# Prepare the same complete-case locked-spec 48h DM 28-day frame
dcc <- lm48_dm_cc
dcc$exp <- dcc$lm_group
# Factor binary covariates so that is.numeric() among covars matches the intended continuous set.
# The prior failure: ck_num_aligned() requires sort(names(which(sapply(df[,covars], is.numeric)))) == sort(num_covar).
# sepsis3/has_hf/has_htn/has_ckd/cardiac_surgery/vent24/vaso24 were stored as numeric 0/1, so they
# were auto-detected as numeric but were omitted from num_covar.
for (v in bin_covars) if (v %in% names(dcc)) dcc[[v]] <- factor(dcc[[v]])
if ("sex" %in% names(dcc) && !is.factor(dcc$sex)) dcc$sex <- factor(dcc$sex)
if ("race_coll" %in% names(dcc) && !is.factor(dcc$race_coll)) dcc$race_coll <- factor(dcc$race_coll)
if ("adm_coll" %in% names(dcc) && !is.factor(dcc$adm_coll)) dcc$adm_coll <- factor(dcc$adm_coll)
if ("year_f" %in% names(dcc) && !is.factor(dcc$year_f)) dcc$year_f <- factor(dcc$year_f)
dcc$lm_group <- factor(dcc$lm_group, levels = c("no_AF", "NOAF", "preexisting"))
num_use <- names(which(sapply(dcc[, covars, drop = FALSE], is.numeric)))
log("aligned num_covar:", paste(num_use, collapse = ", "))

# Always write cox.zph on the identical time-fixed model (TVC PH remains in B07_cox_zph.csv)
f_lm <- as.formula(paste("Surv(lm_time28, status_28) ~ exp +", paste(covars, collapse = " + ")))
fit_lm <- coxph(f_lm, data = dcc, ties = "efron")
zp <- cox.zph(fit_lm)
zph_df <- data.frame(term = rownames(zp$table), chisq = zp$table[, "chisq"],
                     df = zp$table[, "df"], p = zp$table[, "p"],
                     model = "landmark_48h_locked_spec_cc", outcome = "death_28d", subset = "DM",
                     n = nrow(dcc), n_event = as.integer(sum(dcc$status_28)),
                     stringsAsFactors = FALSE)
write_both(zph_df, "B07_cox_zph_landmark48_locked_spec.csv", c("n", "n_event"))
exp_row <- zph_df[zph_df$term == "exp", ]
glob_row <- zph_df[zph_df$term == "GLOBAL", ]
log("cox.zph locked-spec 48h DM 28d: exp p=", if (nrow(exp_row)) exp_row$p else NA,
    "GLOBAL p=", if (nrow(glob_row)) glob_row$p else NA)

strong_note <- data.frame(
  attempted = NA, succeeded = NA, error = "",
  num_covar_used = paste(num_use, collapse = ", "),
  num_covar_intended_continuous = paste(covars_num_intended, collapse = ", "),
  fix = "Converted binary covariates (sepsis3, has_hf, has_htn, has_ckd, cardiac_surgery, vent24, vaso24) to factor so that is.numeric() among covars equals the continuous set; passed that set as num_covar. This satisfies strong:::ck_num_aligned(). The 2026-11-22 add-on call failed because those binaries were still numeric 0/1 in the data frame while num_covar listed only the seven continuous variables.",
  model = "Time-fixed complete-case Cox, locked SAP 48h landmark, diabetes, 28-day death (n should be 12,991). ok_mult_cox does not accept counting-process (tstart,tstop) data, so it cannot be the PH test for the main time-varying models; those remain in B07_cox_zph.csv (survival::cox.zph). op_tdc_tt=FALSE so that strong does not refit with tt(); tbl_pha is only returned when op_tdc_tt=TRUE, therefore PH coefficients for this model are in B07_cox_zph_landmark48_locked_spec.csv (equivalent Schoenfeld test). If tbl_res is written, it is the strong multivariable Cox summary, not a TVC model.",
  stringsAsFactors = FALSE
)

if (!has_strong) {
  strong_note$attempted <- FALSE
  strong_note$succeeded <- FALSE
  strong_note$error <- "strong namespace not available"
  write_both(strong_note, "B07_strong_ok_mult_cox_note.csv")
  log("strong not available; wrote cox.zph equivalent only")
} else {
  strong_note$attempted <- TRUE
  res <- try(strong::ok_mult_cox(
    time = "lm_time28", outcome = "status_28", col_name_gp = "lm_group",
    covars = covars, df = dcc, num_covar = num_use,
    op_tdc_dsp = 1, op_tdc_tt = FALSE,
    p_pha_global = 0.05, p_pha_covar = 0.05,
    mod_label = "landmark48h_locked_spec_DM_28d"
  ), silent = TRUE)
  if (inherits(res, "try-error")) {
    strong_note$succeeded <- FALSE
    strong_note$error <- paste(as.character(res), collapse = " ")
    log("ok_mult_cox still failed:", strong_note$error)
  } else {
    strong_note$succeeded <- TRUE
    strong_note$error <- ""
    # Write tabular pieces only; do not write fit objects
    if (!is.null(res$tbl_res)) {
      tr <- try(as.data.frame(res$tbl_res, stringsAsFactors = FALSE), silent = TRUE)
      if (!inherits(tr, "try-error")) write_both(tr, "B07_strong_ok_mult_cox_tbl_res.csv")
    }
    if (!is.null(res$tbl_pha)) {
      tp <- try(as.data.frame(res$tbl_pha, stringsAsFactors = FALSE), silent = TRUE)
      if (!inherits(tp, "try-error")) write_both(tp, "B07_strong_ok_mult_cox_tbl_pha.csv")
    }
    # tbl_mod may be a gtsummary/gt object; prefer its body dataset if present
    if (!is.null(res$tbl_mod)) {
      tm <- try({
        if (is.data.frame(res$tbl_mod)) as.data.frame(res$tbl_mod)
        else if (!is.null(res$tbl_mod$body$dataset)) as.data.frame(res$tbl_mod$body$dataset)
        else NULL
      }, silent = TRUE)
      if (!inherits(tm, "try-error") && !is.null(tm) && is.data.frame(tm) && nrow(tm) > 0) {
        write_both(tm, "B07_strong_ok_mult_cox_tbl_mod.csv")
      }
    }
    log("ok_mult_cox succeeded; result names:", paste(names(res), collapse = ", "))
  }
  write_both(strong_note, "B07_strong_ok_mult_cox_note.csv")
}

# Confirm official output/ timestamps were not touched (log only)
of <- file.path(proj, "output", "01_flow_counts.csv")
if (file.exists(of)) log("official output/01_flow_counts.csv mtime still", as.character(file.info(of)$mtime))

log("DONE B_addon follow-up 2026-11-22")
