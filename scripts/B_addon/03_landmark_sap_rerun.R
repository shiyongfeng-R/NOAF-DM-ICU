#!/usr/bin/env Rscript
# Rerun 24h and corrected 48h landmarks under SAP spec (alive at landmark via time_28).
# Writes only to output_addon/. No official output/ or results_final/. No individual rows.
# Uses INTERNAL_ONLY _workspace_mids_dat.rds (stay-level analytic frame already built).

suppressPackageStartupMessages({
  library(survival)
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
log_dir <- file.path(proj, "logs")
logf <- file.path(log_dir, "B_addon_landmark_sap.log")
if (file.exists(logf)) file.remove(logf)
log("START landmark SAP-spec rerun 2026-11-22")

n_small <- function(x) { x <- as.numeric(x); !is.na(x) & x >= 1 & x <= 9 }
write_internal <- function(df, name) {
  write.csv(df, file.path(int_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_public <- function(df, name, count_cols = NULL) {
  pub <- df
  if (is.null(count_cols)) {
    count_cols <- names(pub)[grepl("^(n_|n$)|_n$", names(pub)) | names(pub) %in% c("n", "N", "Freq", "n_event")]
  }
  for (cc in count_cols) {
    if (!cc %in% names(df)) next
    x <- as.numeric(df[[cc]])
    pub[[cc]] <- ifelse(is.na(x), "", ifelse(n_small(x), "\u22649",
                     format(as.integer(round(x)), scientific = FALSE, trim = TRUE)))
  }
  write.csv(pub, file.path(out_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_both <- function(df, name, count_cols = NULL) {
  write_internal(df, name); write_public(df, name, count_cols)
}

rds <- file.path(int_dir, "_workspace_mids_dat.rds")
if (!file.exists(rds)) stop("missing INTERNAL_ONLY _workspace_mids_dat.rds")
obj <- readRDS(rds)
dmi0 <- obj$dmi0
# do not print dmi0
log("loaded dmi0 from rds; n not printed")

covars <- c("anchor_age", "sex", "race_coll", "adm_coll", "year_f", "sofa", "sapsii",
            "charlson_ex_dm", "sepsis3", "has_hf", "has_htn", "has_ckd", "cardiac_surgery",
            "vent24", "vaso24", "potassium_min", "mg_min", "gluc_mean")
outcomes <- list(death_28d = c("time_28", "status_28"),
                 death_1y = c("time_1y", "status_1y"),
                 hosp_death = c("time_hosp", "status_hosp"))

subset_dm <- function(d, lab) {
  if (lab == "DM") d[d$has_dm == 1, ] else if (lab == "nonDM") d[d$has_dm == 0, ] else d
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

# require_inhosp: extra still-in-hospital filter used in the earlier add-on
# sap_alive: SAP/locked — alive at landmark by time_28 > horizon (used even for 1-year)
prepare <- function(d, horizon, mode, require_inhosp, oc) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  # SAP alive-at-landmark: time_28 > horizon (locked 48h used this for 28d and 1y)
  d <- d[!is.na(d$time_28) & d$time_28 > horizon, ]
  if (require_inhosp) {
    d <- d[is.na(d$time_hosp) | d$time_hosp > horizon, ]
  }
  if (mode == "exclude_early_noaf") {
    d <- d[!(d$af_group == "NOAF" & !is.na(d$t_af_days) & d$t_af_days <= horizon), ]
    d$lm_group <- as.character(d$af_group)
  } else if (mode == "corrected") {
    d <- d[!(d$af_group == "NOAF" & (is.na(d$t_af_days) | d$t_af_days > horizon)), ]
    d$lm_group <- as.character(d$af_group)
  } else stop("unknown mode")
  d$lm_group <- factor(d$lm_group, levels = c("no_AF", "NOAF", "preexisting"))
  d$lm_time <- pmax(d[[tv]] - horizon, 1 / 24)
  d$lm_status <- d[[sv]]
  d
}

fit_one <- function(d, model, oc, lab) {
  dcc <- d[stats::complete.cases(d[, c("lm_time", "lm_status", "lm_group", covars)]), ]
  if (nrow(dcc) < 50 || sum(dcc$lm_status, na.rm = TRUE) < 10) return(NULL)
  dcc$exp <- dcc$lm_group
  f <- as.formula(paste("Surv(lm_time, lm_status) ~ exp +", paste(covars, collapse = " + ")))
  fit <- try(coxph(f, data = dcc, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  extract_hr(fit, model, oc, lab, n = nrow(dcc), n_event = sum(dcc$lm_status))
}

# ---- HRs: previous extra-filter vs SAP (28d and 1y; DM/nonDM/all) ----
hr_rows <- list()
for (spec in list(
  list(horizon = 1, mode = "exclude_early_noaf", tag = "landmark_1d_exclude_early_noaf"),
  list(horizon = 2, mode = "corrected", tag = "landmark_2d_corrected")
)) {
  for (inhosp in c(TRUE, FALSE)) {
    suf <- if (inhosp) "_addon_inhosp_cc" else "_sap_alive_cc"
    for (lab in c("DM", "nonDM", "all")) {
      for (oc in c("death_28d", "death_1y")) {
        dsub <- subset_dm(dmi0, lab)
        d <- prepare(dsub, spec$horizon, spec$mode, inhosp, oc)
        ex <- fit_one(d, paste0(spec$tag, suf), oc, lab)
        if (!is.null(ex)) {
          ex$hospital_filter <- if (inhosp) "time_hosp > horizon (earlier add-on)" else "none (SAP: alive via time_28)"
          ex$horizon_days <- spec$horizon
          ex$mode <- spec$mode
          hr_rows[[length(hr_rows) + 1]] <- ex
        }
      }
    }
  }
}
hr <- do.call(rbind, hr_rows)
write_both(hr, "B04_landmark_24h_corrected48_sap_rerun.csv", c("n", "n_event", "horizon_days"))
log("wrote SAP vs add-on HR table, rows", nrow(hr))

# ---- Step-by-step counts: DM, 28-day, both landmarks ----
count_steps <- function(horizon, mode, label) {
  d0 <- subset_dm(dmi0, "DM")
  n0 <- nrow(d0)
  ev0 <- sum(d0$status_28 == 1, na.rm = TRUE)
  d_alive <- d0[!is.na(d0$time_28) & d0$time_28 > horizon, ]
  n_alive <- nrow(d_alive)
  still_h <- is.na(d_alive$time_hosp) | d_alive$time_hosp > horizon
  d_inhosp <- d_alive[still_h, ]
  apply_mode <- function(d) {
    if (mode == "exclude_early_noaf") {
      d <- d[!(d$af_group == "NOAF" & !is.na(d$t_af_days) & d$t_af_days <= horizon), ]
      d$lm_group <- as.character(d$af_group)
    } else {
      d <- d[!(d$af_group == "NOAF" & (is.na(d$t_af_days) | d$t_af_days > horizon)), ]
      d$lm_group <- as.character(d$af_group)
    }
    d$lm_group <- factor(d$lm_group, levels = c("no_AF", "NOAF", "preexisting"))
    d$lm_time <- pmax(d$time_28 - horizon, 1 / 24)
    d$lm_status <- d$status_28
    d
  }
  d_sap <- apply_mode(d_alive)
  d_add <- apply_mode(d_inhosp)
  cc <- function(d) d[stats::complete.cases(d[, c("lm_time", "lm_status", "lm_group", covars)]), ]
  d_sap_cc <- cc(d_sap); d_add_cc <- cc(d_add)
  # extra-filter people who were in SAP CC
  extra <- d_alive[!still_h, ]
  extra_m <- apply_mode(extra)
  n_extra_in_sap_cc <- sum(stats::complete.cases(extra_m[, c("lm_time", "lm_status", "lm_group", covars)]))
  gn <- function(d) {
    g <- as.character(d$lm_group)
    c(n_no_AF = sum(g == "no_AF"), n_NOAF = sum(g == "NOAF"), n_pre = sum(g == "preexisting"))
  }
  row <- function(step, spec, desc, d, n_excl, note, verdict = "") {
    if (is.null(d$lm_group)) {
      g <- c(n_no_AF = sum(d$af_group == "no_AF"), n_NOAF = sum(d$af_group == "NOAF"),
             n_pre = sum(d$af_group == "preexisting"))
    } else g <- gn(d)
    data.frame(landmark = label, step = step, spec = spec, description = desc,
               n = as.integer(nrow(d)), n_excluded_this_step = as.integer(n_excl),
               n_event_28 = as.integer(sum(d$status_28 == 1, na.rm = TRUE)),
               n_no_AF = as.integer(g["n_no_AF"]), n_NOAF = as.integer(g["n_NOAF"]),
               n_pre = as.integer(g["n_pre"]), note = note, verdict = verdict,
               stringsAsFactors = FALSE)
  }
  rbind(
    row(1, "common", "Diabetes, three-group, 2008-2019", d0, 0, "Same analytic frame as locked / B00."),
    row(2, "SAP", paste0("Alive at ", horizon * 24, " h (time_28 > ", horizon, ")"),
        d_alive, n0 - n_alive, "SAP / locked definition of being at risk at the landmark. Used for 28-day and 1-year."),
    row(3, "addon_only", "EXTRA: still in hospital (time_hosp > horizon or missing)",
        d_inhosp, n_alive - nrow(d_inhosp),
        "Present in earlier add-on landmark_fit() for 28-day and 1-year. Not in SAP for those outcomes."),
    row(4, "SAP", "Exposure recode/exclusion as specified, then complete-case",
        d_sap_cc, nrow(d_sap) - nrow(d_sap_cc),
        paste0("SAP-spec complete-case n. Of the extra-filter exclusions, ", n_extra_in_sap_cc,
               " were complete-case in the SAP sample and account for the n gap."),
        "CORRECT for 28-day/1-year per SAP (alive via time_28)."),
    row(4, "addon_inhosp", "Same recode/exclusion and complete-case after the extra hospital filter",
        d_add_cc, nrow(d_add) - nrow(d_add_cc),
        "Matches earlier B04_timezero_fixes.csv n for this landmark.",
        "Earlier add-on result; includes a still-in-hospital filter not in the SAP for 28-day/1-year.")
  )
}

recon_new <- rbind(
  count_steps(1, "exclude_early_noaf", "24h_exclude_early_NOAF"),
  count_steps(2, "corrected", "48h_corrected_drop_late_NOAF")
)
# Append to the existing old-style 48h reconciliation (keep prior rows)
old_path <- file.path(int_dir, "B04_landmark48_reconciliation.csv")
if (file.exists(old_path)) {
  old <- read.csv(old_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"landmark" %in% names(old)) old$landmark <- "48h_old_style_recode_late_NOAF_as_noAF"
  # align columns
  cols <- unique(c(names(old), names(recon_new)))
  pad <- function(z) { for (c in setdiff(cols, names(z))) z[[c]] <- ""; z[, cols, drop = FALSE] }
  recon_all <- rbind(pad(old), pad(recon_new))
} else recon_all <- recon_new
write_both(recon_all, "B04_landmark48_reconciliation.csv",
           c("n", "n_excluded_this_step", "n_event_28", "n_no_AF", "n_NOAF", "n_pre"))
log("updated B04_landmark48_reconciliation.csv with 24h and corrected-48h steps")

# Key DM 28-day HRs for the log (aggregates only)
pick <- hr[hr$subset == "DM" & hr$outcome == "death_28d" & hr$term == "NOAF", ]
for (i in seq_len(nrow(pick))) {
  log(pick$model[i], "n=", pick$n[i], "HR", sprintf("%.2f (%.2f-%.2f)", pick$hr[i], pick$lcl[i], pick$ucl[i]))
}

of <- file.path(proj, "output", "01_flow_counts.csv")
if (file.exists(of)) log("official output/01_flow_counts.csv mtime still", as.character(file.info(of)$mtime))
log("DONE landmark SAP-spec rerun")
