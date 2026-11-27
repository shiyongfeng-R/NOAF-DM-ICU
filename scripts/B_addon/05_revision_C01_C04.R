#!/usr/bin/env Rscript
# Methodological revision analyses C01-C04 (DM, 28-day death, TVC Cox, NOAF vs no AF).
# Writes ONLY to output_addon/ (and INTERNAL_ONLY therein).
# Never writes to output/ or results_final/. Never prints individual rows.
# READ_ONLY attach of catalog + mimic.duckdb. No new packages. Seed 20261120.
#
# C01 ECG time-matched confirmation around first_af_time (±1/6/24 h) + HRs
# C02 Pre-AF / time-varying SOFA, vasopressors, invasive ventilation in the main model
# C03 Early (≤48 h) vs late NOAF: baseline + SMD, pre-AF severity aggregates
# C04 Exposure-misclassification grid (Se 0.60-0.95, Sp 0.95-0.99)
#
# Usage (project root):
#   nohup Rscript --vanilla scripts/B_addon/05_revision_C01_C04.R \
#        > logs/B_addon_revision_C.out 2>&1 &

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
logf <- file.path(log_dir, "B_addon_revision_C.log")
if (file.exists(logf)) file.remove(logf)
SEED <- 20261120L
set.seed(SEED)
log("START revision C01-C04; seed", SEED)

n_small <- function(x) { x <- as.numeric(x); !is.na(x) & x >= 1 & x <= 9 }
write_internal <- function(df, name) {
  write.csv(df, file.path(int_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_public <- function(df, name, count_cols = NULL) {
  pub <- df
  if (is.null(count_cols)) {
    count_cols <- names(pub)[grepl("^(n_|n$)|_n$", names(pub)) | names(pub) %in% c("n", "N", "Freq", "n_event")]
  }
  count_cols <- setdiff(count_cols, c("horizon_days", "cut_lo", "cut_hi", "quantile_type",
                                      "se", "sp", "Se", "Sp", "seed", "B", "B_ok"))
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
hr_fmt <- function(hr, lcl, ucl) {
  ifelse(is.na(hr), NA_character_, sprintf("%.2f (%.2f-%.2f)", hr, lcl, ucl))
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
log("analytic frame in memory; in3p defined; n=", sum(in3p), "DM", sum(in3p & dat$has_dm == 1))

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
covars_adm <- c("anchor_age", "sex", "race_coll", "adm_coll", "year_f",
                "charlson_ex_dm", "has_hf", "has_htn", "has_ckd")
form_rhs <- paste(c("exp", covars), collapse = " + ")
imp_vars <- c("sofa", "sapsii", "potassium_min", "mg_min", "gluc_mean")

pool_rubin <- function(coefs, vcovs) {
  m <- length(coefs)
  ubar <- Reduce("+", vcovs) / m
  qbar <- Reduce("+", coefs) / m
  b <- Reduce("+", lapply(coefs, function(q) tcrossprod(q - qbar))) / (m - 1)
  tmat <- ubar + (1 + 1 / m) * b
  se <- sqrt(diag(tmat))
  z <- qbar / se
  data.frame(term = names(qbar), loghr = unname(qbar), se = unname(se),
             hr = unname(exp(qbar)),
             lcl = unname(exp(qbar - 1.96 * se)),
             ucl = unname(exp(qbar + 1.96 * se)),
             p = unname(2 * pnorm(-abs(z))),
             m = m, stringsAsFactors = FALSE)
}

# ---------------- Extra queries on the still-open READ_ONLY connection ----------------
# first_af_time (raw timestamp; t_af_days in dat is floored at 1/48 d)
has_first <- try(DBI::dbExistsTable(con, "analytic") &&
                   "first_af_time" %in% DBI::dbGetQuery(con, "DESCRIBE SELECT * FROM analytic")$column_name,
                 silent = TRUE)
if (!isTRUE(has_first)) {
  # try exposed
  has_first <- try("first_af_time" %in% DBI::dbGetQuery(con, "DESCRIBE SELECT * FROM exposed")$column_name,
                   silent = TRUE)
  src_time <- if (isTRUE(has_first)) "exposed" else NA_character_
} else src_time <- "analytic"
if (is.na(src_time)) {
  log("first_af_time 未确认 in analytic/exposed; will approximate from intime + t_af_days")
  af_times <- data.frame(stay_id = dat$stay_id, first_af_time = as.POSIXct(NA),
                         intime = dat$intime, outtime = dat$outtime, stringsAsFactors = FALSE)
} else {
  af_times <- DBI::dbGetQuery(con, sprintf(
    "SELECT stay_id, first_af_time, intime, outtime FROM %s", src_time))
  log("pulled first_af_time from", src_time, "; n rows not printed")
}
# stay-id lists as in-memory DuckDB temp (not written to official dirs)
stay_dm <- data.frame(stay_id = dat$stay_id[in3p & dat$has_dm == 1])
stay_noaf <- data.frame(stay_id = dat$stay_id[in3p & dat$has_dm == 1 & dat$af_group == "NOAF"])
if (DBI::dbExistsTable(con, "stay_dm")) DBI::dbRemoveTable(con, "stay_dm")
if (DBI::dbExistsTable(con, "stay_noaf")) DBI::dbRemoveTable(con, "stay_noaf")
DBI::dbWriteTable(con, "stay_dm", stay_dm, temporary = TRUE)
DBI::dbWriteTable(con, "stay_noaf", stay_noaf, temporary = TRUE)

# C01: ECG confirmation flags around first_af_time (stay-level flags stay in memory; never printed)
log("C01 SQL: ECG windows around first_af_time")
c01_ok <- TRUE
c01_flags <- try({
  DBI::dbExecute(con, sprintf("
    CREATE OR REPLACE TEMP TABLE c01_flags AS
    SELECT a.stay_id,
      MAX(CASE WHEN e.has_af_any=1 AND e.ecg_time BETWEEN a.first_af_time - INTERVAL 1 HOUR
                 AND a.first_af_time + INTERVAL 1 HOUR THEN 1 ELSE 0 END) AS conf_1h,
      MAX(CASE WHEN e.has_af_any=1 AND e.ecg_time BETWEEN a.first_af_time - INTERVAL 6 HOUR
                 AND a.first_af_time + INTERVAL 6 HOUR THEN 1 ELSE 0 END) AS conf_6h,
      MAX(CASE WHEN e.has_af_any=1 AND e.ecg_time BETWEEN a.first_af_time - INTERVAL 24 HOUR
                 AND a.first_af_time + INTERVAL 24 HOUR THEN 1 ELSE 0 END) AS conf_24h,
      MAX(CASE WHEN e.ecg_time BETWEEN a.first_af_time - INTERVAL 1 HOUR
                 AND a.first_af_time + INTERVAL 1 HOUR THEN 1 ELSE 0 END) AS any_ecg_1h,
      MAX(CASE WHEN e.ecg_time BETWEEN a.first_af_time - INTERVAL 6 HOUR
                 AND a.first_af_time + INTERVAL 6 HOUR THEN 1 ELSE 0 END) AS any_ecg_6h,
      MAX(CASE WHEN e.ecg_time BETWEEN a.first_af_time - INTERVAL 24 HOUR
                 AND a.first_af_time + INTERVAL 24 HOUR THEN 1 ELSE 0 END) AS any_ecg_24h
    FROM %s a
    INNER JOIN stay_noaf s USING (stay_id)
    LEFT JOIN ecg e ON e.subject_id = a.subject_id
    WHERE a.first_af_time IS NOT NULL
    GROUP BY a.stay_id
  ", src_time))
  DBI::dbGetQuery(con, "SELECT * FROM c01_flags")
}, silent = TRUE)
if (inherits(c01_flags, "try-error")) {
  log("C01 ECG window query failed:", as.character(c01_flags)[1])
  c01_ok <- FALSE
  c01_flags <- data.frame(stay_id = integer(), conf_1h = integer(), conf_6h = integer(),
                          conf_24h = integer(), any_ecg_1h = integer(), any_ecg_6h = integer(),
                          any_ecg_24h = integer())
} else log("C01 flags stays", nrow(c01_flags), "(aggregate n only; flags not printed)")

# C02/C03: pre-AF 24h window (and 24-48h) sofa / vent / vaso for DM NOAF; daily TVC for all DM
log("C02/C03 SQL: pre-AF 24h sofa, ventilation, vasopressors")
vaso_ids <- "221906,221289,221662,221653,221749,222315"
# norepinephrine, epinephrine, dopamine, dobutamine, phenylephrine, vasopressin (locked vaso24 set)
preaf <- try({
  DBI::dbExecute(con, sprintf("
    CREATE OR REPLACE TEMP TABLE preaf_sev AS
    SELECT a.stay_id,
      MAX(s.sofa_24hours) FILTER (
        WHERE s.starttime < a.first_af_time
          AND s.endtime > a.first_af_time - INTERVAL 24 HOUR
          AND s.sofa_24hours IS NOT NULL) AS sofa_pre24_max,
      MAX(s.sofa_24hours) FILTER (
        WHERE s.starttime <= a.first_af_time AND s.endtime > a.first_af_time
          AND s.sofa_24hours IS NOT NULL) AS sofa_at_onset,
      MAX(s.sofa_24hours) FILTER (
        WHERE s.starttime < a.first_af_time - INTERVAL 24 HOUR
          AND s.endtime > a.first_af_time - INTERVAL 48 HOUR
          AND s.sofa_24hours IS NOT NULL) AS sofa_pre24_48_max,
      MAX(CASE WHEN v.ventilation_status='InvasiveVent'
                AND v.starttime < a.first_af_time
                AND v.endtime > a.first_af_time - INTERVAL 24 HOUR THEN 1 ELSE 0 END) AS vent_pre24,
      MAX(CASE WHEN v.ventilation_status='InvasiveVent'
                AND v.starttime < a.first_af_time - INTERVAL 24 HOUR
                AND v.endtime > a.first_af_time - INTERVAL 48 HOUR THEN 1 ELSE 0 END) AS vent_pre24_48,
      MAX(CASE WHEN ie.itemid IN (%s)
                AND ie.starttime < a.first_af_time
                AND coalesce(ie.endtime, ie.starttime) > a.first_af_time - INTERVAL 24 HOUR
               THEN 1 ELSE 0 END) AS vaso_pre24,
      MAX(CASE WHEN ie.itemid IN (%s)
                AND ie.starttime < a.first_af_time - INTERVAL 24 HOUR
                AND coalesce(ie.endtime, ie.starttime) > a.first_af_time - INTERVAL 48 HOUR
               THEN 1 ELSE 0 END) AS vaso_pre24_48
    FROM %s a
    INNER JOIN stay_noaf n USING (stay_id)
    LEFT JOIN mimic.mimiciv_derived.sofa s ON s.stay_id=a.stay_id
    LEFT JOIN mimic.mimiciv_derived.ventilation v ON v.stay_id=a.stay_id
    LEFT JOIN wh.mimic_iv.inputevents ie ON ie.stay_id=a.stay_id
    WHERE a.first_af_time IS NOT NULL
    GROUP BY a.stay_id
  ", vaso_ids, vaso_ids, src_time))
  DBI::dbGetQuery(con, "SELECT * FROM preaf_sev")
}, silent = TRUE)
if (inherits(preaf, "try-error")) {
  log("pre-AF severity query failed:", as.character(preaf)[1], "— marked 未确认 where needed")
  preaf <- data.frame(stay_id = integer(), sofa_pre24_max = numeric(), sofa_at_onset = numeric(),
                      sofa_pre24_48_max = numeric(), vent_pre24 = integer(), vent_pre24_48 = integer(),
                      vaso_pre24 = integer(), vaso_pre24_48 = integer())
} else log("pre-AF severity stays", nrow(preaf), "(not printed)")

coalesce_na <- function(x, val = 0) { x[is.na(x)] <- val; x }
log("C02 SQL: daily sofa/vent/vaso for DM 28-day TVC")
daily <- try({
  DBI::dbExecute(con, "
    CREATE OR REPLACE TEMP TABLE days28 AS
    SELECT UNNEST(range(0, 28)) AS day
  ")
  sofa_d <- DBI::dbGetQuery(con, sprintf("
    SELECT a.stay_id,
      CAST(floor(date_diff('second', a.intime, s.starttime)/86400.0) AS INTEGER) AS day,
      MAX(s.sofa_24hours) AS sofa_tv
    FROM mimic.mimiciv_derived.sofa s
    INNER JOIN %s a ON a.stay_id=s.stay_id
    INNER JOIN stay_dm r ON r.stay_id=s.stay_id
    WHERE s.sofa_24hours IS NOT NULL
      AND s.starttime >= a.intime - INTERVAL 1 HOUR
      AND s.starttime < a.intime + INTERVAL 28 DAY
    GROUP BY a.stay_id, day
    HAVING day >= 0 AND day <= 27
  ", src_time))
  vent_d <- DBI::dbGetQuery(con, sprintf("
    SELECT a.stay_id, d.day, 1 AS vent_tv
    FROM mimic.mimiciv_derived.ventilation v
    INNER JOIN %s a ON a.stay_id=v.stay_id
    INNER JOIN stay_dm r ON r.stay_id=v.stay_id
    CROSS JOIN days28 d
    WHERE v.ventilation_status='InvasiveVent'
      AND v.starttime < a.intime + (d.day + 1) * INTERVAL 1 DAY
      AND v.endtime > a.intime + d.day * INTERVAL 1 DAY
      AND v.starttime < a.outtime
    GROUP BY a.stay_id, d.day
  ", src_time))
  vaso_d <- DBI::dbGetQuery(con, sprintf("
    SELECT a.stay_id, d.day, 1 AS vaso_tv
    FROM (
      SELECT ie.stay_id, ie.starttime, coalesce(ie.endtime, ie.starttime) AS endtime
      FROM wh.mimic_iv.inputevents ie
      INNER JOIN stay_dm r ON r.stay_id=ie.stay_id
      WHERE ie.itemid IN (%s)
    ) ie
    INNER JOIN %s a ON a.stay_id=ie.stay_id
    CROSS JOIN days28 d
    WHERE ie.starttime < a.intime + (d.day + 1) * INTERVAL 1 DAY
      AND ie.endtime > a.intime + d.day * INTERVAL 1 DAY
      AND ie.starttime < a.outtime
    GROUP BY a.stay_id, d.day
  ", vaso_ids, src_time))
  # full stay-day grid in R (DM only; never printed)
  grid <- expand.grid(stay_id = stay_dm$stay_id, day = 0:27, stringsAsFactors = FALSE)
  grid <- merge(grid, sofa_d, by = c("stay_id", "day"), all.x = TRUE, sort = FALSE)
  grid <- merge(grid, vent_d, by = c("stay_id", "day"), all.x = TRUE, sort = FALSE)
  grid <- merge(grid, vaso_d, by = c("stay_id", "day"), all.x = TRUE, sort = FALSE)
  grid
}, silent = TRUE)
if (inherits(daily, "try-error")) {
  log("daily severity query failed:", as.character(daily)[1], "— TVC model may be 未确认")
  daily <- NULL
} else {
  log("daily severity rows", nrow(daily), "stays not printed")
  daily$vent_tv <- as.integer(coalesce_na(daily$vent_tv, 0))
  daily$vaso_tv <- as.integer(coalesce_na(daily$vaso_tv, 0))
}

# LOCF sofa by stay (after pull; do not print)
if (!is.null(daily) && nrow(daily)) {
  daily <- daily[order(daily$stay_id, daily$day), ]
  # simple LOCF within stay
  last <- NA_real_; last_id <- daily$stay_id[1] - 1L
  so <- daily$sofa_tv
  sid <- daily$stay_id
  for (i in seq_along(so)) {
    if (sid[i] != last_id) { last <- NA_real_; last_id <- sid[i] }
    if (!is.na(so[i])) last <- so[i] else so[i] <- last
  }
  daily$sofa_tv_locf <- so
}

# close DB — no official writes were made
try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
log("disconnected DuckDB (READ_ONLY session ended)")

# merge stay-level extras onto a working copy (never printed / never written as rows)
dat$af_group <- factor(dat$af_group, levels = c("no_AF", "NOAF", "preexisting", "grayzone", "indeterminate"))
work <- dat[in3p, ]
work <- merge(work, af_times[, c("stay_id", "first_af_time")], by = "stay_id", all.x = TRUE, sort = FALSE)
if (nrow(c01_flags)) {
  work <- merge(work, c01_flags, by = "stay_id", all.x = TRUE, sort = FALSE)
} else {
  work$conf_1h <- work$conf_6h <- work$conf_24h <- NA_integer_
  work$any_ecg_1h <- work$any_ecg_6h <- work$any_ecg_24h <- NA_integer_
}
if (nrow(preaf)) {
  work <- merge(work, preaf, by = "stay_id", all.x = TRUE, sort = FALSE)
} else {
  work$sofa_pre24_max <- work$sofa_at_onset <- work$sofa_pre24_48_max <- NA_real_
  work$vent_pre24 <- work$vent_pre24_48 <- work$vaso_pre24 <- work$vaso_pre24_48 <- NA_integer_
}
# for non-NOAF, pre-AF window is undefined: use first-day sofa/vent24/vaso24
is_noaf <- as.character(work$af_group) == "NOAF"
work$sofa_pre24 <- ifelse(is_noaf, work$sofa_at_onset, work$sofa)
work$sofa_pre24[is_noaf & is.na(work$sofa_pre24)] <- work$sofa_pre24_max[is_noaf & is.na(work$sofa_pre24)]
# fallback to first-day sofa if still missing
work$sofa_pre24[is.na(work$sofa_pre24)] <- work$sofa[is.na(work$sofa_pre24)]
work$vent_pre24_use <- ifelse(is_noaf, work$vent_pre24, work$vent24)
work$vaso_pre24_use <- ifelse(is_noaf, work$vaso_pre24, work$vaso24)
work$vent_pre24_use[is.na(work$vent_pre24_use)] <- work$vent24[is.na(work$vent_pre24_use)]
work$vaso_pre24_use[is.na(work$vaso_pre24_use)] <- work$vaso24[is.na(work$vaso_pre24_use)]

# ---------------- MICE (same specification as locked / B_addon) ----------------
log("mice m=10 maxit=8 seed", SEED, "on in3p covariates only")
dmi0 <- work
pred_df <- dmi0[, c("has_dm", "af_group", covars)]
for (vv in c("sepsis3", "has_hf", "has_htn", "has_ckd", "cardiac_surgery", "vent24", "vaso24")) {
  pred_df[[vv]] <- factor(pred_df[[vv]])
}
do_mice <- TRUE
ini <- try(mice::mice(pred_df, maxit = 0, print = FALSE), silent = TRUE)
if (inherits(ini, "try-error")) {
  log("mice setup failed"); do_mice <- FALSE; imp <- NULL
} else {
  meth <- ini$method
  for (v in names(meth)) meth[v] <- ""
  for (v in imp_vars) if (v %in% names(meth)) meth[v] <- "pmm"
  imp <- try(mice::mice(pred_df, m = 10, maxit = 8, method = meth, seed = SEED, print = FALSE), silent = TRUE)
  if (inherits(imp, "try-error")) { log("mice failed"); do_mice <- FALSE; imp <- NULL }
  else log("mice completed")
}
fill_imp <- function(k) {
  dfill <- dmi0
  if (isTRUE(do_mice)) {
    comp <- mice::complete(imp, k)
    for (v in imp_vars) dfill[[v]] <- comp[[v]]
    # if sofa_pre24 fell back to first-day sofa and that was missing, refresh from imputed sofa
    miss_pre <- is.na(dfill$sofa_pre24) | (as.character(dfill$af_group) != "NOAF" & is.na(dmi0$sofa))
    # for non-NOAF, sofa_pre24 tracks first-day sofa
    dfill$sofa_pre24[as.character(dfill$af_group) != "NOAF"] <- dfill$sofa[as.character(dfill$af_group) != "NOAF"]
    still <- is.na(dfill$sofa_pre24)
    dfill$sofa_pre24[still] <- dfill$sofa[still]
  }
  dfill
}
subset_dm <- function(d, lab = "DM") {
  if (lab == "DM") d[d$has_dm == 1, ] else if (lab == "nonDM") d[d$has_dm == 0, ] else d
}

fit_tvc_rhs <- function(d, rhs) {
  d <- drop_noaf_missing_onset(d, "tvc")
  long <- restore_factors(make_tvc(d, "time_28", "status_28"))
  need <- unique(c("tstart", "tstop", "event", "exp", intersect(all.vars(as.formula(paste("~", rhs))), names(long))))
  long <- long[stats::complete.cases(long[, need]), ]
  if (nrow(long) < 50 || sum(long$event) < 10) return(NULL)
  f <- as.formula(paste("Surv(tstart, tstop, event) ~", rhs))
  fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}
fit_mi_rhs <- function(d0, rhs, model, drop_ids = NULL) {
  coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
  for (k in 1:10) {
    dfill <- subset_dm(fill_imp(k), "DM")
    if (!is.null(drop_ids)) dfill <- dfill[!dfill$stay_id %in% drop_ids, ]
    fit <- fit_tvc_rhs(dfill, rhs)
    if (is.null(fit)) next
    cf <- coef(fit); vc <- vcov(fit)
    if (!"expNOAF" %in% names(cf)) next
    keep <- grepl("^exp", names(cf))
    coefs[[length(coefs) + 1]] <- cf[keep]
    vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
    ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill$status_28, na.rm = TRUE))
  }
  if (length(coefs) < 3) return(NULL)
  po <- pool_rubin(coefs, vcovs)
  po <- po[po$term == "expNOAF", , drop = FALSE]
  if (!nrow(po)) return(NULL)
  po$term <- "NOAF"
  po$model <- model
  po$outcome <- "death_28d"
  po$subset <- "DM"
  po$contrast_vs <- "no_AF"
  po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
  po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
  po$m_ok <- length(coefs)
  po
}

# =====================================================================
# C01 ECG time-matched confirmation
# =====================================================================
log("C01: ECG time-matched confirmation and HRs")
ddm <- subset_dm(dmi0, "DM")
n_noaf <- sum(as.character(ddm$af_group) == "NOAF")
n_noaf_time <- sum(as.character(ddm$af_group) == "NOAF" & !is.na(ddm$first_af_time))
n_noaf_miss_time <- n_noaf - n_noaf_time
c01_counts <- list()
c01_hrs <- list()
for (w in c(1, 6, 24)) {
  conf_col <- paste0("conf_", w, "h")
  any_col <- paste0("any_ecg_", w, "h")
  noaf <- ddm[as.character(ddm$af_group) == "NOAF", ]
  n_conf <- sum(noaf[[conf_col]] == 1, na.rm = TRUE)
  n_any <- sum(noaf[[any_col]] == 1, na.rm = TRUE)
  n_miss_flag <- sum(as.character(ddm$af_group) == "NOAF" & is.na(ddm[[conf_col]]))
  c01_counts[[length(c01_counts) + 1]] <- data.frame(
    window = paste0("\u00b1", w, " h"),
    window_hours = w,
    n_NOAF_DM = n_noaf,
    n_NOAF_with_first_af_time = n_noaf_time,
    n_NOAF_missing_first_af_time = n_noaf_miss_time,
    n_confirmed_AF_ECG = n_conf,
    n_any_ECG_in_window = n_any,
    n_flag_missing = n_miss_flag,
    prop_confirmed_among_NOAF = if (n_noaf > 0) n_conf / n_noaf else NA_real_,
    prop_confirmed_among_NOAF_with_ECG_in_window = if (n_any > 0) n_conf / n_any else NA_real_,
    note = if (!c01_ok) "未确认: ECG window query failed"
           else if (n_noaf_miss_time > 0)
             "Confirmation uses first_af_time (chart 220048). Stays missing first_af_time cannot be confirmed in a time window (counted in denominator of prop_confirmed_among_NOAF, not in n_confirmed)."
           else "Anchor = first Heart Rhythm (itemid 220048) AF/AFL charttime (first_af_time). Confirmed = at least one machine-read ECG (same regex as locked pipeline) with AF/AFL within the window.",
    stringsAsFactors = FALSE
  )
  # HR: exclude unconfirmed NOAF (do not recode them as no_AF)
  if (c01_ok && n_conf >= 10) {
    drop_ids <- noaf$stay_id[!(noaf[[conf_col]] == 1)]
    po <- try(fit_mi_rhs(dmi0, form_rhs, paste0("tvc_mice_ecg_confirmed_pm", w, "h"),
                         drop_ids = drop_ids), silent = TRUE)
    if (!inherits(po, "try-error") && !is.null(po)) {
      po$window <- paste0("\u00b1", w, " h")
      po$n_NOAF_confirmed <- n_conf
      po$n_NOAF_excluded_unconfirmed <- n_noaf - n_conf
      po$note <- "Time-varying Cox, MI m=10, same covariates as the locked main model. Unconfirmed NOAF stays excluded (not recoded as no AF). Pre-existing AF retained as the third exposure level. Contrast = confirmed NOAF vs no AF."
      c01_hrs[[length(c01_hrs) + 1]] <- po
      log("C01 window", w, "h confirmed", n_conf, "HR", po$hr_ci)
    } else log("C01 HR failed for window", w)
  } else log("C01 HR skipped for window", w, "confirmed n=", n_conf)
}
write_both(do.call(rbind, c01_counts), "C01_ecg_time_match_counts.csv",
           c("n_NOAF_DM", "n_NOAF_with_first_af_time", "n_NOAF_missing_first_af_time",
             "n_confirmed_AF_ECG", "n_any_ECG_in_window", "n_flag_missing"))
if (length(c01_hrs)) {
  write_both(do.call(rbind, c01_hrs), "C01_ecg_time_match_hr.csv",
             c("n", "n_event", "n_NOAF_confirmed", "n_NOAF_excluded_unconfirmed", "m", "m_ok"))
} else {
  write_both(data.frame(item = "C01_hr", status = "未确认",
                        reason = "ECG-confirmed HR could not be estimated (query failed or too few confirmed events)",
                        stringsAsFactors = FALSE), "C01_ecg_time_match_hr.csv")
}
log("C01 done")

# =====================================================================
# C02 pre-AF / time-varying severity in the main model
# =====================================================================
log("C02: pre-AF and time-varying severity models")
# Definitions (for README and a lookup table)
c02_def <- data.frame(
  variable = c("sofa (main model)", "vent24 (main model)", "vaso24 (main model)",
               "sofa_pre24", "vent_pre24", "vaso_pre24",
               "sofa_tv", "vent_tv", "vaso_tv"),
  definition = c(
    "First-day SOFA from mimiciv_derived.first_day_sofa (locked main model).",
    "Any invasive ventilation (mimiciv_derived.ventilation, ventilation_status='InvasiveVent') overlapping the first 24 h after ICU intime (locked).",
    "Any vasopressor start in the first 24 h after ICU intime, inputevents itemids 221906 (norepinephrine), 221289 (epinephrine), 221662 (dopamine), 221653 (dobutamine), 221749 (phenylephrine), 222315 (vasopressin) (locked vaso24).",
    "For NOAF: sofa_24hours from mimiciv_derived.sofa in the hourly interval containing first_af_time (fallback: max sofa_24hours in the 24 h ending at first_af_time; then first-day SOFA). For no AF / pre-existing AF: first-day SOFA (pre-AF window is undefined).",
    "For NOAF: any InvasiveVent interval overlapping (first_af_time-24 h, first_af_time]. For no AF / pre-existing AF: vent24.",
    "For NOAF: any of the locked vaso itemids with infusion overlapping (first_af_time-24 h, first_af_time]. For no AF / pre-existing AF: vaso24.",
    "Daily (calendar day from ICU intime, days 0-27): max sofa_24hours among mimiciv_derived.sofa intervals overlapping that day; last observation carried forward after the last ICU value, including after ICU discharge. If no hourly SOFA is ever available, first-day (imputed) SOFA is used.",
    "1 if an InvasiveVent interval overlaps that day and starts before ICU outtime; 0 after ICU discharge or if none.",
    "1 if a locked vaso itemid infusion overlaps that day and starts before ICU outtime; 0 after ICU discharge or if none."
  ),
  source_table = c("mimic.mimiciv_derived.first_day_sofa",
                   "mimic.mimiciv_derived.ventilation",
                   "wh.mimic_iv.inputevents (itemids as above)",
                   "mimic.mimiciv_derived.sofa (sofa_24hours, starttime, endtime)",
                   "mimic.mimiciv_derived.ventilation",
                   "wh.mimic_iv.inputevents",
                   "mimic.mimiciv_derived.sofa",
                   "mimic.mimiciv_derived.ventilation",
                   "wh.mimic_iv.inputevents"),
  stringsAsFactors = FALSE
)
write_both(c02_def, "C02_variable_definitions.csv")

c02_rows <- list()
# A. Locked main model, refit (should be ~1.67)
po_main <- try(fit_mi_rhs(dmi0, form_rhs, "tvc_mice_main_refit"), silent = TRUE)
if (!inherits(po_main, "try-error") && !is.null(po_main)) {
  po_main$comparator <- "main_model_target_1.67"
  po_main$note <- "Refit of the locked main TVC MI model in the diabetes 28-day cohort (same covariates)."
  c02_rows[[length(c02_rows) + 1]] <- po_main
  log("C02 main refit", po_main$hr_ci)
}
# B. Admission-only (target 2.48) — refit for the same session
rhs_adm <- paste(c("exp", covars_adm), collapse = " + ")
po_adm <- try(fit_mi_rhs(dmi0, rhs_adm, "tvc_mice_admission_covars_refit"), silent = TRUE)
if (!inherits(po_adm, "try-error") && !is.null(po_adm)) {
  po_adm$comparator <- "admission_only_target_2.48"
  po_adm$note <- "Refit of the add-on admission-only TVC (demographics, comorbidities, admission type, year; no first-day severity/treatment/labs)."
  c02_rows[[length(c02_rows) + 1]] <- po_adm
  log("C02 admission-only refit", po_adm$hr_ci)
}
# C. Replace sofa, vent24, vaso24 with pre-AF 24 h values
rhs_pre <- paste(c("exp", "anchor_age", "sex", "race_coll", "adm_coll", "year_f",
                   "sofa_pre24", "sapsii", "charlson_ex_dm", "sepsis3", "has_hf", "has_htn", "has_ckd",
                   "cardiac_surgery", "vent_pre24_use", "vaso_pre24_use",
                   "potassium_min", "mg_min", "gluc_mean"), collapse = " + ")
po_pre <- try(fit_mi_rhs(dmi0, rhs_pre, "tvc_mice_preaf24_severity"), silent = TRUE)
if (!inherits(po_pre, "try-error") && !is.null(po_pre)) {
  po_pre$comparator <- "vs_main_1.67_and_admission_2.48"
  po_pre$note <- "Locked covariate set, but sofa/vent24/vaso24 replaced by pre-AF 24 h values for NOAF and by first-day values for no AF and pre-existing AF."
  c02_rows[[length(c02_rows) + 1]] <- po_pre
  log("C02 pre-AF 24h", po_pre$hr_ci)
} else log("C02 pre-AF 24h model 未确认")

# D. Time-varying sofa, vent, vaso (drop first-day sofa/vent24/vaso24 to avoid duplication)
fit_tvc_daily <- function(d, daily) {
  d <- drop_noaf_missing_onset(d, "tvc_daily")
  long <- restore_factors(make_tvc(d, "time_28", "status_28"))
  # attach first-day sofa as day-0 fallback
  fb <- d[, c("stay_id", "sofa")]
  names(fb)[2] <- "sofa_day0"
  daily2 <- merge(daily, fb, by = "stay_id", all.x = TRUE, sort = FALSE)
  daily2$sofa_use <- daily2$sofa_tv_locf
  daily2$sofa_use[is.na(daily2$sofa_use)] <- daily2$sofa_day0[is.na(daily2$sofa_use)]
  daily2$vent_use <- as.integer(coalesce_na(daily2$vent_tv, 0))
  daily2$vaso_use <- as.integer(coalesce_na(daily2$vaso_tv, 0))
  mer <- merge(long, daily2[, c("stay_id", "day", "sofa_use", "vent_use", "vaso_use")],
               by = "stay_id", all.x = TRUE, sort = FALSE)
  mer <- mer[!is.na(mer$day) & (mer$day + 1) > mer$tstart & mer$day < mer$tstop, ]
  mer$tstart2 <- pmax(mer$tstart, mer$day)
  mer$tstop2 <- pmin(mer$tstop, mer$day + 1)
  mer <- mer[mer$tstop2 > mer$tstart2, ]
  # event only on the last sub-interval of the original row
  mer$event2 <- 0L
  mer <- mer[order(mer$stay_id, mer$tstart, mer$tstart2), ]
  # last row within each original interval (stay_id, tstart, tstop)
  key <- paste(mer$stay_id, mer$tstart, mer$tstop)
  last <- !duplicated(key, fromLast = TRUE)
  mer$event2[last] <- mer$event[last]
  mer$tstart <- mer$tstart2; mer$tstop <- mer$tstop2; mer$event <- mer$event2
  rhs <- paste(c("exp", "anchor_age", "sex", "race_coll", "adm_coll", "year_f",
                 "sofa_use", "sapsii", "charlson_ex_dm", "sepsis3", "has_hf", "has_htn", "has_ckd",
                 "cardiac_surgery", "vent_use", "vaso_use",
                 "potassium_min", "mg_min", "gluc_mean"), collapse = " + ")
  need <- unique(c("tstart", "tstop", "event", "exp", "sofa_use", "vent_use", "vaso_use",
                   "anchor_age", "sex", "race_coll", "adm_coll", "year_f", "sapsii",
                   "charlson_ex_dm", "sepsis3", "has_hf", "has_htn", "has_ckd",
                   "cardiac_surgery", "potassium_min", "mg_min", "gluc_mean"))
  mer <- restore_factors(mer)
  mer <- mer[stats::complete.cases(mer[, need]), ]
  if (nrow(mer) < 50 || sum(mer$event) < 10) return(NULL)
  f <- as.formula(paste("Surv(tstart, tstop, event) ~", rhs))
  fit <- try(coxph(f, data = mer, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}
if (!is.null(daily) && nrow(daily)) {
  coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
  for (k in 1:10) {
    dfill <- subset_dm(fill_imp(k), "DM")
    fit <- try(fit_tvc_daily(dfill, daily), silent = TRUE)
    if (inherits(fit, "try-error") || is.null(fit)) next
    cf <- coef(fit); vc <- vcov(fit)
    if (!"expNOAF" %in% names(cf)) next
    keep <- grepl("^exp", names(cf))
    coefs[[length(coefs) + 1]] <- cf[keep]
    vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
    ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill$status_28, na.rm = TRUE))
  }
  if (length(coefs) >= 3) {
    po <- pool_rubin(coefs, vcovs)
    po <- po[po$term == "expNOAF", , drop = FALSE]
    po$term <- "NOAF"
    po$model <- "tvc_mice_severity_tv_daily"
    po$outcome <- "death_28d"; po$subset <- "DM"; po$contrast_vs <- "no_AF"
    po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
    po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
    po$m_ok <- length(coefs)
    po$comparator <- "vs_main_1.67_and_admission_2.48"
    po$note <- "Locked covariate set, but first-day sofa/vent24/vaso24 replaced by daily time-varying sofa_24hours, invasive ventilation, and vasopressors (LOCF sofa after last ICU value; vent/vaso 0 after ICU discharge)."
    c02_rows[[length(c02_rows) + 1]] <- po
    log("C02 daily TVC", po$hr_ci)
  } else log("C02 daily TVC 未确认 (too few MI fits)")
} else log("C02 daily TVC 未确认 (no daily severity table)")

if (length(c02_rows)) {
  # pad columns
  nms <- unique(unlist(lapply(c02_rows, names)))
  c02_rows <- lapply(c02_rows, function(z) { for (m in setdiff(nms, names(z))) z[[m]] <- NA; z[, nms, drop = FALSE] })
  c02 <- do.call(rbind, c02_rows)
  write_both(c02, "C02_preaf_severity_hr.csv", c("n", "n_event", "m", "m_ok"))
} else {
  write_both(data.frame(item = "C02", status = "未确认", stringsAsFactors = FALSE),
             "C02_preaf_severity_hr.csv")
}
log("C02 done")

# =====================================================================
# C03 early (≤48 h) vs late NOAF
# =====================================================================
log("C03: early vs late NOAF baseline and pre-AF trajectory")
noaf <- subset_dm(dmi0, "DM")
noaf <- noaf[as.character(noaf$af_group) == "NOAF", ]
noaf <- drop_noaf_missing_onset(noaf, "C03")
noaf$early <- as.integer(!is.na(noaf$t_af_days) & noaf$t_af_days <= 2)
n_early <- sum(noaf$early == 1); n_late <- sum(noaf$early == 0)
log("C03 n early", n_early, "late", n_late)

smd_cont <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  x1 <- x1[!is.na(x1)]; x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  (mean(x1) - mean(x0)) / sqrt((stats::var(x1) + stats::var(x0)) / 2)
}
smd_bin <- function(x, g) {
  x <- as.numeric(x)
  p1 <- mean(x[g == 1], na.rm = TRUE); p0 <- mean(x[g == 0], na.rm = TRUE)
  if (is.na(p1) || is.na(p0)) return(NA_real_)
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
}
row_cont <- function(label, x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  data.frame(variable = label, type = "continuous",
             n_early = sum(!is.na(x1)), n_late = sum(!is.na(x0)),
             early_mean = mean(x1, na.rm = TRUE), early_sd = stats::sd(x1, na.rm = TRUE),
             early_median = stats::median(x1, na.rm = TRUE),
             early_p25 = unname(stats::quantile(x1, 0.25, na.rm = TRUE)),
             early_p75 = unname(stats::quantile(x1, 0.75, na.rm = TRUE)),
             late_mean = mean(x0, na.rm = TRUE), late_sd = stats::sd(x0, na.rm = TRUE),
             late_median = stats::median(x0, na.rm = TRUE),
             late_p25 = unname(stats::quantile(x0, 0.25, na.rm = TRUE)),
             late_p75 = unname(stats::quantile(x0, 0.75, na.rm = TRUE)),
             n_early_yes = NA_integer_, n_late_yes = NA_integer_,
             early_pct = NA_real_, late_pct = NA_real_,
             smd = smd_cont(x, g), stringsAsFactors = FALSE)
}
row_bin <- function(label, x, g) {
  x <- as.numeric(x)
  data.frame(variable = label, type = "binary",
             n_early = sum(g == 1 & !is.na(x)), n_late = sum(g == 0 & !is.na(x)),
             early_mean = NA_real_, early_sd = NA_real_, early_median = NA_real_,
             early_p25 = NA_real_, early_p75 = NA_real_,
             late_mean = NA_real_, late_sd = NA_real_, late_median = NA_real_,
             late_p25 = NA_real_, late_p75 = NA_real_,
             n_early_yes = sum(x[g == 1] == 1, na.rm = TRUE),
             n_late_yes = sum(x[g == 0] == 1, na.rm = TRUE),
             early_pct = 100 * mean(x[g == 1] == 1, na.rm = TRUE),
             late_pct = 100 * mean(x[g == 0] == 1, na.rm = TRUE),
             smd = smd_bin(x, g), stringsAsFactors = FALSE)
}
g <- noaf$early
tab <- rbind(
  row_cont("anchor_age", noaf$anchor_age, g),
  row_bin("male", as.integer(noaf$sex == "M"), g),
  row_bin("race_White", as.integer(noaf$race_coll == "White"), g),
  row_bin("admission_Emergency", as.integer(noaf$adm_coll == "Emergency"), g),
  row_cont("sofa_first_day", noaf$sofa, g),
  row_cont("sapsii", noaf$sapsii, g),
  row_cont("charlson_ex_dm", noaf$charlson_ex_dm, g),
  row_bin("sepsis3", noaf$sepsis3, g),
  row_bin("has_hf", noaf$has_hf, g),
  row_bin("has_htn", noaf$has_htn, g),
  row_bin("has_ckd", noaf$has_ckd, g),
  row_bin("cardiac_surgery", noaf$cardiac_surgery, g),
  row_bin("vent24", noaf$vent24, g),
  row_bin("vaso24", noaf$vaso24, g),
  row_cont("potassium_min", noaf$potassium_min, g),
  row_cont("mg_min", noaf$mg_min, g),
  row_cont("gluc_mean", noaf$gluc_mean, g),
  row_cont("t_af_days", noaf$t_af_days, g),
  row_cont("icu_los_days", noaf$icu_los_days, g),
  row_bin("death_28d", noaf$status_28, g)
)
tab$n_early_group <- n_early
tab$n_late_group <- n_late
tab$definition_early <- "t_af_days <= 2 (48 h from ICU intime)"
write_both(tab, "C03_early_late_baseline.csv",
           c("n_early", "n_late", "n_early_yes", "n_late_yes", "n_early_group", "n_late_group"))

# Trajectory aggregates (no individual paths)
traj_block <- function(label, sofa, vent, vaso, g, note) {
  rbind(
    { r <- row_cont(paste0(label, "_sofa"), sofa, g); r$metric_note <- note; r },
    { r <- row_bin(paste0(label, "_vent"), vent, g); r$metric_note <- note; r },
    { r <- row_bin(paste0(label, "_vaso"), vaso, g); r$metric_note <- note; r }
  )
}
traj <- rbind(
  traj_block("first24h", noaf$sofa, noaf$vent24, noaf$vaso24, g,
             "ICU first 24 h (locked first-day SOFA, vent24, vaso24)."),
  traj_block("preAF_24to0h", noaf$sofa_pre24_max, noaf$vent_pre24, noaf$vaso_pre24, g,
             "24 h window ending at first_af_time: max sofa_24hours; any invasive vent; any vasopressor. 未确认 if first_af_time or the source interval is missing."),
  traj_block("preAF_at_onset", noaf$sofa_at_onset, noaf$vent_pre24, noaf$vaso_pre24, g,
             "sofa_24hours in the hourly SOFA interval containing first_af_time; vent/vaso as in the 24 h window ending at onset."),
  traj_block("preAF_48to24h", noaf$sofa_pre24_48_max, noaf$vent_pre24_48, noaf$vaso_pre24_48, g,
             "Window 48-24 h before first_af_time. For early NOAF (onset ≤48 h) this window often overlaps pre-ICU time and is frequently missing (report n; do not impute).")
)
traj$n_early_group <- n_early
traj$n_late_group <- n_late
write_both(traj, "C03_preaf_severity_trajectory.csv",
           c("n_early", "n_late", "n_early_yes", "n_late_yes", "n_early_group", "n_late_group"))
log("C03 done")

# =====================================================================
# C04 misclassification grid
# =====================================================================
log("C04: Se/Sp misclassification grid, complete-case TVC, B=200, seed", SEED)
# Sample: DM, 28-day, complete-case on locked covars; preexisting held fixed
cc_vars <- c("time_28", "status_28", "af_group", covars)
dcc <- subset_dm(dmi0, "DM")
dcc <- dcc[stats::complete.cases(dcc[, cc_vars]), ]
dcc <- drop_noaf_missing_onset(dcc, "C04")
dcc$af_group <- factor(as.character(dcc$af_group), levels = c("no_AF", "NOAF", "preexisting"))
n_cc <- nrow(dcc)
n_noaf_cc <- sum(dcc$af_group == "NOAF")
n_noaf_grp <- sum(dcc$af_group %in% c("no_AF", "NOAF"))
p_hat <- if (n_noaf_grp > 0) n_noaf_cc / n_noaf_grp else NA_real_
log("C04 CC n=", n_cc, "NOAF", n_noaf_cc, "p_hat among no_AF+NOAF", sprintf("%.4f", p_hat))
t_af_pool <- dcc$t_af_days[dcc$af_group == "NOAF" & !is.na(dcc$t_af_days)]

fit_one_cc <- function(d) {
  long <- restore_factors(make_tvc(d, "time_28", "status_28"))
  long <- long[stats::complete.cases(long[, c("tstart", "tstop", "event", "exp", covars)]), ]
  if (nrow(long) < 50 || sum(long$event) < 10) return(NA_real_)
  f <- as.formula(paste("Surv(tstart, tstop, event) ~", form_rhs))
  fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  cf <- coef(fit)
  if (!"expNOAF" %in% names(cf)) return(NA_real_)
  unname(cf["expNOAF"])
}

# observed CC HR once (for the method note)
b_obs <- fit_one_cc(dcc)
hr_obs <- if (is.finite(b_obs)) exp(b_obs) else NA_real_
log("C04 observed CC HR NOAF vs no_AF", hr_obs)

se_grid <- seq(0.60, 0.95, by = 0.05)
sp_grid <- seq(0.95, 0.99, by = 0.01)
B_mc <- 200L
set.seed(SEED)
# Pre-draw uniform randoms so each cell is reproducible from the same seed stream
# We use a cell-specific subseed: SEED + 1000*i_se + 10*i_sp

c04_rows <- list()
for (i in seq_along(se_grid)) {
  for (j in seq_along(sp_grid)) {
    Se <- se_grid[i]; Sp <- sp_grid[j]
    PPV <- Se * p_hat / (Se * p_hat + (1 - Sp) * (1 - p_hat))
    NPV <- Sp * (1 - p_hat) / ((1 - Se) * p_hat + Sp * (1 - p_hat))
    subseed <- as.integer(SEED + 1000L * i + 10L * j)
    set.seed(subseed)
    hs <- rep(NA_real_, B_mc)
    id_noaf <- which(dcc$af_group == "NOAF")
    id_no   <- which(dcc$af_group == "no_AF")
    for (b in seq_len(B_mc)) {
      x <- dcc
      # FP: classified NOAF, truly no AF, with prob 1-PPV
      if (length(id_noaf)) {
        fp <- stats::runif(length(id_noaf)) > PPV
        if (any(fp)) {
          x$af_group[id_noaf[fp]] <- "no_AF"
          x$t_af_days[id_noaf[fp]] <- NA_real_
        }
      }
      # FN: classified no AF, truly NOAF, with prob 1-NPV; onset sampled from observed NOAF t_af
      if (length(id_no) && length(t_af_pool)) {
        fn <- stats::runif(length(id_no)) > NPV
        if (any(fn)) {
          x$af_group[id_no[fn]] <- "NOAF"
          ta <- sample(t_af_pool, sum(fn), replace = TRUE)
          # keep onset strictly before follow-up end
          tmax <- pmax(x$time_28[id_no[fn]] - 1 / 48, 1 / 48)
          ta <- pmin(ta, tmax)
          x$t_af_days[id_no[fn]] <- ta
        }
      }
      x$af_group <- factor(as.character(x$af_group), levels = c("no_AF", "NOAF", "preexisting"))
      hs[b] <- fit_one_cc(x)
    }
    ok <- is.finite(hs)
    c04_rows[[length(c04_rows) + 1]] <- data.frame(
      Se = Se, Sp = Sp, PPV = PPV, NPV = NPV,
      p_hat = p_hat, n = n_cc, n_NOAF = n_noaf_cc, n_noAF_plus_NOAF = n_noaf_grp,
      B = B_mc, B_ok = sum(ok), seed_base = SEED, seed_cell = subseed,
      hr_observed_cc = hr_obs,
      hr_median = if (any(ok)) exp(stats::median(hs[ok])) else NA_real_,
      hr_mean = if (any(ok)) exp(mean(hs[ok])) else NA_real_,
      hr_lcl = if (any(ok)) exp(stats::quantile(hs[ok], 0.025)) else NA_real_,
      hr_ucl = if (any(ok)) exp(stats::quantile(hs[ok], 0.975)) else NA_real_,
      stringsAsFactors = FALSE
    )
    log("C04 Se", Se, "Sp", Sp, "median HR",
        if (any(ok)) sprintf("%.2f", exp(stats::median(hs[ok]))) else "NA",
        "B_ok", sum(ok))
  }
}
c04 <- do.call(rbind, c04_rows)
write_both(c04, "C04_misclassification_grid.csv",
           c("n", "n_NOAF", "n_noAF_plus_NOAF", "B", "B_ok"))
c04_method <- data.frame(
  item = "C04",
  method = paste(
    "Record-level Monte Carlo correction for non-differential misclassification of NOAF vs no AF",
    "in the diabetes, 28-day, complete-case time-varying Cox model (same covariates as the locked main model).",
    "Pre-existing AF is held fixed and is not reclassified.",
    "Prevalence p_hat is the observed NOAF share among no-AF + NOAF in the complete-case sample.",
    "PPV = Se*p / (Se*p + (1-Sp)*(1-p)); NPV = Sp*(1-p) / ((1-Se)*p + Sp*(1-p)).",
    "Each iteration: classified NOAF are recoded to no AF with probability 1-PPV (onset cleared);",
    "classified no AF are recoded to NOAF with probability 1-NPV, with onset time sampled with replacement",
    "from the observed NOAF t_af_days in the same sample, capped below follow-up end.",
    "The sampled-onset rule for false negatives is an assumption (true missed-onset times are 未确认).",
    "B=200 iterations per (Se, Sp) cell; report median and 2.5th-97.5th percentiles of the HR (exp of log HR).",
    "Base seed 20261120; cell seed = 20261120 + 1000*i_Se + 10*j_Sp (i,j 1-based indices on the grids).",
    "Se grid 0.60-0.95 by 0.05; Sp grid 0.95-0.99 by 0.01.",
    "This extends the earlier PPV-only QBA (B05) and is not a multiple-imputation pooled HR."
  ),
  stringsAsFactors = FALSE
)
write_both(c04_method, "C04_method_note.csv")
log("C04 done, cells", nrow(c04))

# session
sess <- data.frame(
  item = c("date", "timezone", "R", "survival", "mice", "duckdb", "seed", "cohort",
           "C01", "C02", "C03", "C04"),
  value = c(as.character(Sys.time()), "Pacific/Apia (UTC+13)",
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(utils::packageVersion("survival")),
            as.character(utils::packageVersion("mice")),
            as.character(utils::packageVersion("duckdb")),
            as.character(SEED),
            "DM, 2008-2019, three-group in3p; outcome 28-day death; contrast NOAF vs no AF, TVC Cox",
            "ECG confirmation around first_af_time",
            "pre-AF 24h and daily time-varying sofa/vent/vaso",
            "early vs late baseline + SMD and pre-AF aggregates",
            "Se/Sp grid record-level MC, B=200"),
  stringsAsFactors = FALSE
)
write_both(sess, "C00_revision_session_info.csv")

of <- file.path(proj, "output", "01_flow_counts.csv")
if (file.exists(of)) log("official output/01_flow_counts.csv mtime still", as.character(file.info(of)$mtime))
log("DONE revision C01-C04")
