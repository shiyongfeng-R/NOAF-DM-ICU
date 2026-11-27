#!/usr/bin/env Rscript
# Cosmetic fixes to add-on aggregate CSVs (2026-11-27).
# Writes only under output_addon/ (and INTERNAL_ONLY therein).
# Does not write to output/ or results_final/. Does not print individual rows.
# Does not install packages. Does not edit manuscripts.
#
# 1. B11: highest-tertile cut_lo was wrongly set to the lower cut; set it to the
#    upper tertile cut (glucose CV 0.2446..., mean glucose 186 mg/dL).
# 2. SAP rerun public file: horizon_days (1 or 2) is not a count and must not be
#    masked as <=9. Restore from INTERNAL_ONLY; other count-masking rules unchanged.
# 3. B04_timezero_fixes.csv: mark landmark rows that used the still-in-hospital
#    filter as superseded; do not delete them.

suppressPackageStartupMessages({
  # base + utils only; no new packages
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
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(log_dir, "B_addon_fix_20261127.log")
if (file.exists(logf)) file.remove(logf)
log <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = logf, append = TRUE)
}
log("START add-on output typo/masking/annotation fixes")

if (!dir.exists(out_dir) || !dir.exists(int_dir)) stop("output_addon paths missing")

n_small <- function(x) {
  x <- as.numeric(x)
  !is.na(x) & x >= 1 & x <= 9
}
write_internal <- function(df, name) {
  write.csv(df, file.path(int_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_public <- function(df, name, count_cols = NULL) {
  pub <- df
  if (is.null(count_cols)) {
    count_cols <- names(pub)[grepl("^(n_|n$)|_n$", names(pub)) | names(pub) %in% c("n", "N", "Freq", "n_event")]
  }
  # never treat these as counts
  count_cols <- setdiff(count_cols, c("horizon_days", "cut_lo", "cut_hi", "quantile_type"))
  for (cc in count_cols) {
    if (!cc %in% names(df)) next
    x <- as.numeric(df[[cc]])
    pub[[cc]] <- ifelse(is.na(x), "",
                        ifelse(n_small(x), "\u22649",
                               format(as.integer(round(x)), scientific = FALSE, trim = TRUE)))
  }
  write.csv(pub, file.path(out_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_both <- function(df, name, count_cols = NULL) {
  write_internal(df, name)
  write_public(df, name, count_cols)
}

# ---------------- 1. B11 highest-tertile cut_lo ----------------
# Prefer INTERNAL (unrounded source values); fall back to public if needed.
b11_path <- file.path(int_dir, "B11_tableS8_tertile_cutpoints.csv")
if (!file.exists(b11_path)) b11_path <- file.path(out_dir, "B11_tableS8_tertile_cutpoints.csv")
b11 <- read.csv(b11_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
# For each variable, upper cut = middle tertile cut_hi (same as locked quantile 2/3)
fix_hi <- function(df, varname, high_tertile, mid_tertile) {
  i_hi <- which(df$variable == varname & df$tertile == high_tertile)
  i_mid <- which(df$variable == varname & df$tertile == mid_tertile)
  if (length(i_hi) != 1L || length(i_mid) != 1L) {
    stop("could not uniquely locate tertile rows for ", varname)
  }
  old <- df$cut_lo[i_hi]
  new <- df$cut_hi[i_mid]
  df$cut_lo[i_hi] <- new
  log("B11", varname, high_tertile, "cut_lo", old, "->", new, "; cut_hi remains empty (open upper bound); interval already correct")
  df
}
b11 <- fix_hi(b11, "gluc_cv", "CVt3_high", "CVt2")
b11 <- fix_hi(b11, "gluc_mean", "meant3_high", "meant2")
write_both(b11, "B11_tableS8_tertile_cutpoints.csv",
           c("n_with_variable", "n_complete_case", "n_event", "n_missing_variable_in_DM_in3p"))

# ---------------- 2. Restore horizon_days in the public SAP rerun file ----------------
sap_int <- file.path(int_dir, "B04_landmark_24h_corrected48_sap_rerun.csv")
if (!file.exists(sap_int)) stop("missing INTERNAL SAP rerun file")
sap <- read.csv(sap_int, stringsAsFactors = FALSE, check.names = FALSE)
# horizon_days must be the integer 1 or 2, not a masked count
if (!"horizon_days" %in% names(sap)) stop("horizon_days column missing in INTERNAL file")
# coerce to integer (INTERNAL has the true values)
sap$horizon_days <- as.integer(sap$horizon_days)
if (any(!sap$horizon_days %in% c(1L, 2L))) {
  log("WARNING unexpected horizon_days values after restore")
}
write_both(sap, "B04_landmark_24h_corrected48_sap_rerun.csv",
           c("n", "n_event"))  # deliberately omit horizon_days
log("SAP rerun: restored horizon_days to 1/2 in public file; n/n_event masking unchanged")

# ---------------- 3. Annotate superseded landmark rows in B04_timezero_fixes ----------------
tz_int <- file.path(int_dir, "B04_timezero_fixes.csv")
if (!file.exists(tz_int)) tz_int <- file.path(out_dir, "B04_timezero_fixes.csv")
tz <- read.csv(tz_int, stringsAsFactors = FALSE, check.names = FALSE)
lm_models <- c("landmark_1d_exclude_early_noaf_cc",
               "landmark_2d_corrected_cc",
               "landmark_2d_old_style_cc")
is_lm <- tz$model %in% lm_models
tz$status <- ifelse(is_lm, "superseded", "current")
tz$status_note <- ifelse(is_lm,
  "Superseded: this landmark used a still-in-hospital filter (time_hosp > horizon) that is not in the SAP for 28-day/1-year outcomes. Use B04_landmark_24h_corrected48_sap_rerun.csv (SAP: alive via time_28) and B04_landmark48_reconciliation.csv. Rows retained for the record.",
  "")
write_both(tz, "B04_timezero_fixes.csv", c("n", "n_event"))
log("timezero_fixes: marked", sum(is_lm), "of", nrow(tz), "rows as superseded (landmark models); other rows current; none deleted")

# Confirm official output was not touched
of <- file.path(proj, "output", "01_flow_counts.csv")
if (file.exists(of)) log("official output/01_flow_counts.csv mtime still", as.character(file.info(of)$mtime))
log("DONE add-on output typo/masking/annotation fixes")
