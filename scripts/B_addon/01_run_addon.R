#!/usr/bin/env Rscript
# CVD add-on analyses (items 1-10). Aggregates and figures only.
# Writes ONLY to output_addon/ (and INTERNAL_ONLY therein).
# Never writes to output/ or results_final/. Never prints individual rows.
# Does not install packages. Seed 20261120 (SAP).
#
# Packages used:
#   survival, mice, ggplot2, duckdb, DBI  -- locked pipeline
#   cmprsk, tidycmprsk (if available)    -- competing-risk CIF
#   strong (v1.0.9, optional at runtime) -- fhcv/fpv formatting; ok_mult_cox
#       PHA table for time-fixed landmark models only (does not support
#       counting-process / Simon-Makuch / RERI / QBA)
#   mimicR 1.0 is a PostgreSQL MIMIC connector (not used; locked pipeline
#       uses duckdb READ_ONLY attach). mimicR430 has km_plot / stratum_model
#       for baseline-time Cox/KM, not time-varying NOAF.
#
# Usage (from anywhere):
#   nohup Rscript --vanilla scripts/B_addon/01_run_addon.R \
#        > logs/B_addon.out 2>&1 &

suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(survival); library(ggplot2)
})

# ---------------- paths: source helpers then immediately redirect I/O ----------------
# Paths come from environment variables / auto-detection (see scripts/00_helpers.R and README.md).
source(if (file.exists("scripts/00_helpers.R")) "scripts/00_helpers.R" else {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (!length(f)) stop("Set NOAF_REPO and run from the repository root (README.md).")
  file.path(dirname(dirname(dirname(normalizePath(f[1])))), "scripts/00_helpers.R")
}, local = FALSE)
# 00_helpers.R sets out_dir to output/; override BEFORE any write or log to main log
addon_dir <- file.path(proj, "output_addon")
out_dir <- addon_dir
int_dir <- file.path(addon_dir, "INTERNAL_ONLY")
fig_dir <- file.path(addon_dir, "figures")
log_dir <- file.path(proj, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(log_dir, "B_addon.log")
if (file.exists(logf)) file.remove(logf)
SEED <- 20261120L
set.seed(SEED)
log("START B_addon CVD supplementary analyses; seed", SEED)

# Optional strong formatters (package requires network licence check; do not abort if missing)
has_strong <- isTRUE(requireNamespace("strong", quietly = TRUE))
log("strong available:", has_strong, if (has_strong) as.character(packageVersion("strong")) else "")
fh2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_character_, sprintf("%.2f", x))
}
fp2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), NA_character_,
         ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}
if (has_strong) {
  try({
    fh2 <<- function(x) {
      x <- suppressWarnings(as.numeric(x))
      ifelse(is.na(x), NA_character_, as.character(strong::fhcv(x)))
    }
    fp2 <<- function(x) {
      x <- suppressWarnings(as.numeric(x))
      ifelse(is.na(x), NA_character_, {
        z <- strong::fpv(x, equal_sign = FALSE)
        as.character(z)
      })
    }
    log("using strong::fhcv / strong::fpv for table formatting")
  }, silent = TRUE)
}

# ---------------- writers (mask n<=9 in public files) ----------------
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
    if (cc %in% c("n", "N", "n_total", "n_event", "n_atrisk")) row_denom_small <- row_denom_small | sm
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

hr_fmt <- function(hr, lcl, ucl) {
  ifelse(is.na(hr), NA_character_, paste0(fh2(hr), " (", fh2(lcl), "-", fh2(ucl), ")"))
}

bind_fill <- function(...) {
  xs <- list(...)
  if (length(xs) == 1 && is.list(xs[[1]]) && !is.data.frame(xs[[1]])) xs <- xs[[1]]
  xs <- xs[!vapply(xs, is.null, logical(1))]
  xs <- xs[vapply(xs, function(z) is.data.frame(z) && nrow(z) > 0, logical(1))]
  if (!length(xs)) return(NULL)
  nms <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(z) {
    miss <- setdiff(nms, names(z))
    for (m in miss) z[[m]] <- NA
    z[, nms, drop = FALSE]
  })
  do.call(rbind, xs)
}

# ---------------- Item 10 can be written from warehouse metadata (no rows) ----------------
item10 <- data.frame(
  item = "eICU_CRD_version",
  version = "2.0",
  source = "CoredataWarehouse/metadata",
  evidence = paste(
    "metadata/datasets.csv: eicu, PhysioNet eICU-CRD, version 2.0;",
    "metadata/raw_paths.tsv: PhysioNet eICU-CRD 2.0; do not modify;",
    "metadata/attach_paths.tsv: eicu duckdb is attached READ_ONLY; set NOAF_EICU_DUCKDB for a local file;",
    "analysis uses catalog views wh.eicu_parquet (31-table parquet set);",
    "patient table n=200859 stays, 208 hospitals (official eICU-CRD 2.0 anchors);",
    "SAP_v1.md also names eICU-CRD 2.0.",
    sep = " "
  ),
  stringsAsFactors = FALSE
)
write_both(item10, "B10_eicu_version.csv")
log("wrote B10_eicu_version.csv : 2.0")

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
# sanity aggregates only
san <- data.frame(
  n_dat = nrow(dat),
  n_in3p = sum(in3p),
  n_in3p_dm = sum(in3p & dat$has_dm == 1),
  n_NOAF_0819 = sum(in3p & dat$af_group == "NOAF"),
  n_NOAF_dm = sum(in3p & dat$af_group == "NOAF" & dat$has_dm == 1),
  n_pre = sum(in3p & dat$af_group == "preexisting"),
  n_noAF = sum(in3p & dat$af_group == "no_AF"),
  stringsAsFactors = FALSE
)
write_both(san, "B00_prep_sanity_counts.csv")
log("prep sanity three-group 2008-19 n=", san$n_in3p, "DM", san$n_in3p_dm, "NOAF", san$n_NOAF_0819)

# functions copied from locked 01_run_main.R (defined after Flow in the original)
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
form_rhs <- paste(c("exp", covars), collapse = " + ")
covars_num <- c("anchor_age", "sofa", "sapsii", "charlson_ex_dm", "potassium_min", "mg_min", "gluc_mean")
# Admission-only (no first-24h severity / treatment / labs)
covars_adm <- c("anchor_age", "sex", "race_coll", "adm_coll", "year_f",
                "charlson_ex_dm", "has_hf", "has_htn", "has_ckd")
# Admission-only plus comorbidities that may partly reflect early ICU course
covars_adm2 <- c(covars_adm, "cardiac_surgery", "sepsis3")

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

restore_factors <- function(d) {
  if ("sex" %in% names(d)) d$sex <- factor(d$sex, levels = c("F", "M"))
  if ("race_coll" %in% names(d)) d$race_coll <- factor(d$race_coll, levels = c("White", "Black", "Hispanic", "Asian", "Other or Unknown"))
  if ("adm_coll" %in% names(d)) d$adm_coll <- factor(d$adm_coll, levels = c("Emergency", "Elective", "Other"))
  if ("year_f" %in% names(d)) d$year_f <- factor(d$year_f)
  if ("exp" %in% names(d) && !is.factor(d$exp)) d$exp <- factor(d$exp)
  d
}

subset_dm <- function(d, lab) {
  if (lab == "DM") d[d$has_dm == 1, , drop = FALSE]
  else if (lab == "nonDM") d[d$has_dm == 0, , drop = FALSE]
  else d
}

# Rubin pool of a named coefficient, or of a contrast (diff of two coefs)
pool_rubin <- function(coefs, vcovs, contrast = NULL) {
  # coefs, vcovs: lists of length m; contrast = c(name_a, name_b) for b_a - b_b, or a single name
  tn <- if (is.null(contrast)) names(coefs[[1]]) else contrast
  out <- list()
  if (is.null(contrast) || length(contrast) == 1) {
    terms <- if (is.null(contrast)) tn else contrast
    for (tm in terms) {
      qs <- sapply(coefs, function(b) unname(b[tm]))
      if (any(is.na(qs))) next
      ubar <- mean(sapply(seq_along(vcovs), function(i) vcovs[[i]][tm, tm]))
      bvar <- stats::var(qs)
      m <- length(qs)
      Tvar <- ubar + (1 + 1 / m) * bvar
      est <- mean(qs); se <- sqrt(Tvar)
      out[[length(out) + 1]] <- data.frame(term = sub("^exp", "", tm),
        loghr = est, se = se, hr = exp(est),
        lcl = exp(est - 1.96 * se), ucl = exp(est + 1.96 * se),
        p = 2 * pnorm(-abs(est / se)), m = m, stringsAsFactors = FALSE)
    }
  } else {
    a <- contrast[1]; b <- contrast[2]
    qs <- sapply(coefs, function(cf) unname(cf[a] - cf[b]))
    ubar <- mean(sapply(seq_along(vcovs), function(i) {
      v <- vcovs[[i]]
      v[a, a] + v[b, b] - 2 * v[a, b]
    }))
    bvar <- stats::var(qs)
    m <- length(qs)
    Tvar <- ubar + (1 + 1 / m) * bvar
    est <- mean(qs); se <- sqrt(Tvar)
    out[[1]] <- data.frame(term = "NOAF_vs_preexisting",
      loghr = est, se = se, hr = exp(est),
      lcl = exp(est - 1.96 * se), ucl = exp(est + 1.96 * se),
      p = 2 * pnorm(-abs(est / se)), m = m, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

fit_tvc_one <- function(d, time_var, status_var, rhs = form_rhs) {
  d <- drop_noaf_missing_onset(d)
  if (nrow(d) < 50 || sum(d[[status_var]], na.rm = TRUE) < 10) return(NULL)
  long <- restore_factors(make_tvc(d, time_var, status_var))
  f <- as.formula(paste("Surv(tstart, tstop, event) ~", rhs))
  fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

# ---------------- Item 9: exact covariate list from locked code ----------------
item9 <- data.frame(
  role = c(rep("exposure", 1), rep("covariate_in_form_rhs", length(covars)),
           "imputed_variables", "not_in_imputation_model", "mice_settings"),
  name = c(
    "exp (time-varying: no_AF / NOAF / preexisting; ties='efron')",
    covars,
    paste(c("sofa", "sapsii", "potassium_min", "mg_min", "gluc_mean"), collapse = ", "),
    "outcome indicator and Nelson-Aalen cumulative hazard are NOT in pred_df",
    "mice PMM, m=10, maxit=8, seed=20261120; predictors = has_dm + af_group + covars"
  ),
  in_SAP_s6 = c(
    "yes",
    "age yes", "sex yes", "race yes (collapsed)", "admission type yes",
    "year group yes", "first-day SOFA yes", "SAPS II yes",
    "Charlson excluding diabetes (Quan, index admission ICD)",
    "Sepsis-3 derived table (not ICD-based)",
    "HF ICD yes", "HTN ICD yes", "CKD ICD yes",
    "cardiac surgery (CSURG service or CABG/valve ICD) yes",
    "invasive vent in first 24h yes", "vasoactive in first 24h yes",
    "potassium min first 24h yes", "magnesium min first 24h yes", "mean glucose first 24h yes",
    "yes (only these five are imputed)",
    "SAP said impute baseline covariates only, not exposure or outcome; Nelson-Aalen was not implemented",
    "matches SAP m and seed"
  ),
  stringsAsFactors = FALSE
)
write_both(item9, "B09_main_model_covariates.csv")
log("wrote B09_main_model_covariates.csv; race, adm type, Charlson ex-DM, year_f ARE in the locked model")

# ---------------- mice (same specification as 01_run_main.R) ----------------
log("mice m=10 maxit=8 seed", SEED, "on in3p covariates only")
imp_vars <- c("sofa", "sapsii", "potassium_min", "mg_min", "gluc_mean")
dmi0 <- dat[in3p, c("stay_id", "has_dm", "af_group", "t_af_days",
                    "time_hosp", "status_hosp", "time_28", "status_28", "time_1y", "status_1y",
                    "has_ecg_hosp", "af_hosp_ecg", "prior_af_icd", "af_pre", "early_af_2h",
                    "icu_los_days", "died_in_icu", "days_to_discharge", "days_to_death",
                    covars)]
# do not print dmi0
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
# stay-level objects for reruns stay on the Mac only
if (isTRUE(do_mice)) {
  try(saveRDS(list(dmi0 = dmi0, imp = imp, in3p_n = san$n_in3p, seed = SEED),
              file = file.path(int_dir, "_workspace_mids_dat.rds"), compress = "xz"), silent = TRUE)
}

fill_imp <- function(k) {
  dfill <- dmi0
  if (isTRUE(do_mice)) {
    comp <- mice::complete(imp, k)
    for (v in imp_vars) dfill[[v]] <- comp[[v]]
  }
  dfill
}

outcomes <- list(
  death_28d = c("time_28", "status_28"),
  hosp_death = c("time_hosp", "status_hosp"),
  death_1y = c("time_1y", "status_1y")
)
dm_labs <- c("DM", "nonDM", "all")

# ---------------- Item 1: NOAF vs pre-existing (MI pooled, with CI from vcov) ----------------
log("Item 1: NOAF vs preexisting, MI-pooled contrast")
item1_rows <- list()
for (oc in names(outcomes)) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  for (lab in dm_labs) {
    coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
    for (k in 1:10) {
      dfill <- subset_dm(fill_imp(k), lab)
      fit <- fit_tvc_one(dfill, tv, sv)
      if (is.null(fit)) next
      cf <- coef(fit); vc <- vcov(fit)
      keep <- grepl("^exp", names(cf))
      coefs[[length(coefs) + 1]] <- cf[keep]
      vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
      ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill[[sv]], na.rm = TRUE))
    }
    if (length(coefs) < 3) { log("too few MI fits", oc, lab); next }
    vs_no <- pool_rubin(coefs, vcovs, contrast = NULL)
    vs_pre <- pool_rubin(coefs, vcovs, contrast = c("expNOAF", "exppreexisting"))
    for (i in seq_len(nrow(vs_no))) {
      item1_rows[[length(item1_rows) + 1]] <- data.frame(
        model = "tvc_mice_m10", outcome = oc, subset = lab,
        term = vs_no$term[i], contrast_vs = "no_AF",
        hr = vs_no$hr[i], lcl = vs_no$lcl[i], ucl = vs_no$ucl[i], p = vs_no$p[i],
        hr_ci = hr_fmt(vs_no$hr[i], vs_no$lcl[i], vs_no$ucl[i]),
        n = round(mean(ns)), n_event = round(mean(evs)),
        stringsAsFactors = FALSE)
    }
    item1_rows[[length(item1_rows) + 1]] <- data.frame(
      model = "tvc_mice_m10", outcome = oc, subset = lab,
      term = "NOAF", contrast_vs = "preexisting",
      hr = vs_pre$hr, lcl = vs_pre$lcl, ucl = vs_pre$ucl, p = vs_pre$p,
      hr_ci = hr_fmt(vs_pre$hr, vs_pre$lcl, vs_pre$ucl),
      n = round(mean(ns)), n_event = round(mean(evs)),
      stringsAsFactors = FALSE)
  }
}
item1 <- bind_fill(item1_rows)
# compare A1 point estimates
a1_map <- data.frame(
  outcome = c("death_28d", "death_28d", "death_28d", "hosp_death", "hosp_death", "hosp_death",
              "death_1y", "death_1y", "death_1y"),
  subset = c("DM", "nonDM", "all", "DM", "nonDM", "all", "DM", "nonDM", "all"),
  a1_point = c(1.1496, 1.1941, 1.1872, 1.034, 1.0674, 1.0682, 0.9711, 1.1062, 1.063),
  stringsAsFactors = FALSE
)
item1$a1_point <- NA_real_
for (i in seq_len(nrow(item1))) {
  if (item1$contrast_vs[i] == "preexisting" && item1$term[i] == "NOAF") {
    j <- which(a1_map$outcome == item1$outcome[i] & a1_map$subset == item1$subset[i])
    if (length(j) == 1) item1$a1_point[i] <- a1_map$a1_point[j]
  }
}
item1$a1_minus_verified <- item1$hr - item1$a1_point
write_both(item1, "B01_noaf_vs_preexisting.csv", c("n", "n_event"))
log("Item 1 done, rows", nrow(item1))

# ---------------- Item 2: Knol RERI / AP / S (CC + delta + bootstrap; optional MI) ----------------
log("Item 2: Knol-recoded RERI (ref = DM and no NOAF)")
# 4-level: A=NOAF, B=no diabetes; A0B0 = DM & no NOAF
reri_from_hrs <- function(h10, h01, h11) {
  # h10 = HR_A1B0 (NOAF, DM); h01 = HR_A0B1 (no NOAF, nonDM); h11 = HR_A1B1 (NOAF, nonDM)
  reri <- h11 - h10 - h01 + 1
  ap <- reri / h11
  s <- (h11 - 1) / ((h10 - 1) + (h01 - 1))
  c(HR_NOAF_DM = h10, HR_noNOAF_nonDM = h01, HR_NOAF_nonDM = h11, RERI = reri, AP = ap, S = s)
}
delta_reri <- function(fit) {
  # grp4 levels A0B0, A1B0, A0B1, A1B1; coefs are log HR vs A0B0
  cf <- coef(fit); V <- vcov(fit)
  nms <- c("grp4A1B0", "grp4A0B1", "grp4A1B1")
  if (!all(nms %in% names(cf))) return(NULL)
  b10 <- unname(cf["grp4A1B0"]); b01 <- unname(cf["grp4A0B1"]); b11 <- unname(cf["grp4A1B1"])
  h10 <- exp(b10); h01 <- exp(b01); h11 <- exp(b11)
  est <- reri_from_hrs(h10, h01, h11)
  # RERI = e^b11 - e^b10 - e^b01 + 1
  g_reri <- c(-h10, -h01, h11)
  # AP = RERI/h11 = 1 - e^{b10-b11} - e^{b01-b11} + e^{-b11}
  # d AP / db10 = -e^{b10-b11}; d AP/db01 = -e^{b01-b11}; d AP/db11 = e^{b10-b11}+e^{b01-b11}-e^{-b11}
  g_ap <- c(-exp(b10 - b11), -exp(b01 - b11), exp(b10 - b11) + exp(b01 - b11) - exp(-b11))
  # S = (h11-1)/((h10-1)+(h01-1)); more awkward; use numerical gradient if needed
  den <- (h10 - 1) + (h01 - 1)
  # dS/db10 = -h10*(h11-1)/den^2 ; dS/db01 = -h01*(h11-1)/den^2 ; dS/db11 = h11/den
  g_s <- c(-h10 * (h11 - 1) / den^2, -h01 * (h11 - 1) / den^2, h11 / den)
  VV <- V[nms, nms, drop = FALSE]
  se_reri <- sqrt(as.numeric(t(g_reri) %*% VV %*% g_reri))
  se_ap <- sqrt(as.numeric(t(g_ap) %*% VV %*% g_ap))
  se_s <- sqrt(as.numeric(t(g_s) %*% VV %*% g_s))
  data.frame(
    HR_NOAF_DM = h10, HR_noNOAF_nonDM = h01, HR_NOAF_nonDM = h11,
    RERI = est["RERI"], RERI_lcl = est["RERI"] - 1.96 * se_reri, RERI_ucl = est["RERI"] + 1.96 * se_reri,
    AP = est["AP"], AP_lcl = est["AP"] - 1.96 * se_ap, AP_ucl = est["AP"] + 1.96 * se_ap,
    S = est["S"], S_lcl = est["S"] - 1.96 * se_s, S_ucl = est["S"] + 1.96 * se_s,
    method = "delta", stringsAsFactors = FALSE
  )
}

make_grp4_tvc <- function(d, time_var, status_var) {
  # exclude preexisting; 2x2 of NOAF x diabetes, then recode
  d <- d[d$af_group %in% c("no_AF", "NOAF"), ]
  d <- drop_noaf_missing_onset(d, "reri")
  long <- make_tvc(d, time_var, status_var)
  long$grp4 <- factor(
    paste0(ifelse(as.character(long$exp) == "NOAF", "A1", "A0"),
           ifelse(long$has_dm == 0, "B1", "B0")),
    levels = c("A0B0", "A1B0", "A0B1", "A1B1")
  )
  restore_factors(long)
}

fit_reri_cc <- function(d, time_var, status_var, rhs_cov = paste(covars, collapse = " + ")) {
  long <- make_grp4_tvc(d, time_var, status_var)
  long <- long[stats::complete.cases(long[, c("tstart", "tstop", "event", "grp4", covars)]), ]
  if (nrow(long) < 50) return(list(fit = NULL, n = 0, n_event = 0, n_stay = 0))
  f <- as.formula(paste("Surv(tstart, tstop, event) ~ grp4 +", rhs_cov))
  fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(fit = NULL, n = 0, n_event = 0, n_stay = 0))
  list(fit = fit, n = length(unique(long$stay_id)), n_event = sum(long$event), n_stay = length(unique(long$stay_id)), long = long)
}

reri_boot <- function(d, time_var, status_var, B = 500, seed = SEED) {
  set.seed(seed)
  ids <- unique(d$stay_id)
  nst <- length(ids)
  row_by_id <- split(seq_len(nrow(d)), d$stay_id)
  mat <- matrix(NA_real_, nrow = B, ncol = 6)
  colnames(mat) <- c("HR_NOAF_DM", "HR_noNOAF_nonDM", "HR_NOAF_nonDM", "RERI", "AP", "S")
  for (b in seq_len(B)) {
    take <- sample(ids, nst, replace = TRUE)
    idx <- unlist(row_by_id[as.character(take)], use.names = FALSE)
    db <- d[idx, , drop = FALSE]
    # new stay ids for resampled copies (aggregates only; never printed)
    times <- vapply(row_by_id[as.character(take)], length, integer(1))
    db$stay_id <- rep(seq_along(take), times = times)
    ft <- fit_reri_cc(db, time_var, status_var)
    if (is.null(ft$fit)) next
    dl <- delta_reri(ft$fit)
    if (is.null(dl)) next
    mat[b, ] <- c(dl$HR_NOAF_DM, dl$HR_noNOAF_nonDM, dl$HR_NOAF_nonDM, dl$RERI, dl$AP, dl$S)
    if (b %% 50 == 0) log("  bootstrap iter", b, "/", B)
  }
  ok <- stats::complete.cases(mat)
  list(mat = mat, n_ok = sum(ok),
       est = colMeans(mat[ok, , drop = FALSE]),
       lcl = apply(mat[ok, , drop = FALSE], 2, quantile, 0.025, na.rm = TRUE),
       ucl = apply(mat[ok, , drop = FALSE], 2, quantile, 0.975, na.rm = TRUE))
}

# bind data frames that may have different columns
bind_fill <- function(...) {
  xs <- list(...)
  # also accept a list
  if (length(xs) == 1 && is.list(xs[[1]]) && !is.data.frame(xs[[1]])) xs <- xs[[1]]
  xs <- xs[!vapply(xs, is.null, logical(1))]
  xs <- xs[vapply(xs, function(z) is.data.frame(z) && nrow(z) > 0, logical(1))]
  if (!length(xs)) return(NULL)
  nms <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(z) {
    miss <- setdiff(nms, names(z))
    for (m in miss) z[[m]] <- NA
    z[, nms, drop = FALSE]
  })
  do.call(rbind, xs)
}

# Complete-case interaction cohort (same rule as locked RERI)
int_base <- dat[in3p & dat$af_group %in% c("no_AF", "NOAF"), ]
item2_rows <- list()
a2_map <- list(
  death_28d = list(RERI = 0.1549, AP = 0.0639, S = 1.122, n = 38039),
  hosp_death = list(RERI = -0.0008, AP = -0.0004, S = 0.9993, n = 38039),
  death_1y = list(RERI = 0.2081, AP = 0.1232, S = 1.4323, n = 38039)
)
B_boot <- c(death_28d = 500, hosp_death = 300, death_1y = 300)
for (oc in names(outcomes)) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  dcc <- int_base[stats::complete.cases(int_base[, c(tv, sv, "af_group", "has_dm", covars)]), ]
  dcc <- drop_noaf_missing_onset(dcc, paste("RERI CC", oc))
  ft <- fit_reri_cc(dcc, tv, sv)
  if (is.null(ft$fit)) { log("RERI CC fit failed", oc); next }
  dl <- delta_reri(ft$fit)
  log("RERI CC", oc, "n_stay", ft$n, "delta RERI", round(dl$RERI, 4))
  item2_rows[[length(item2_rows) + 1]] <- cbind(
    outcome = oc, model = "tvc_interaction_cc", n = ft$n, n_event = ft$n_event,
    dl, a2_RERI = a2_map[[oc]]$RERI, a2_AP = a2_map[[oc]]$AP, a2_S = a2_map[[oc]]$S,
    stringsAsFactors = FALSE)
  # bootstrap (stay-level)
  log("RERI bootstrap", oc, "B", B_boot[[oc]])
  bt <- try(reri_boot(dcc, tv, sv, B = B_boot[[oc]], seed = SEED), silent = TRUE)
  if (!inherits(bt, "try-error")) {
    item2_rows[[length(item2_rows) + 1]] <- data.frame(
      outcome = oc, model = "tvc_interaction_cc_bootstrap", n = ft$n, n_event = ft$n_event,
      HR_NOAF_DM = bt$est["HR_NOAF_DM"], HR_noNOAF_nonDM = bt$est["HR_noNOAF_nonDM"], HR_NOAF_nonDM = bt$est["HR_NOAF_nonDM"],
      RERI = bt$est["RERI"], RERI_lcl = bt$lcl["RERI"], RERI_ucl = bt$ucl["RERI"],
      AP = bt$est["AP"], AP_lcl = bt$lcl["AP"], AP_ucl = bt$ucl["AP"],
      S = bt$est["S"], S_lcl = bt$lcl["S"], S_ucl = bt$ucl["S"],
      method = paste0("bootstrap_B", B_boot[[oc]], "_ok", bt$n_ok),
      a2_RERI = a2_map[[oc]]$RERI, a2_AP = a2_map[[oc]]$AP, a2_S = a2_map[[oc]]$S,
      stringsAsFactors = FALSE)
  } else log("bootstrap failed", oc)
}
item2 <- bind_fill(item2_rows)
write_both(item2, "B02_knol_reri.csv", c("n", "n_event"))
log("Item 2 done")

# ---------------- Item 4 helpers used also below: pre-existing split ----------------
dmi0$pre_documented <- as.integer(dmi0$prior_af_icd == 1 | dmi0$af_pre == 1)
dmi0$af_at_admission_only <- as.integer(as.character(dmi0$af_group) == "preexisting" &
                                          dmi0$prior_af_icd == 0 & dmi0$af_pre == 0)
# counts
split_n <- aggregate(stay_id ~ has_dm + I(as.character(af_group)) + pre_documented + af_at_admission_only,
                     data = dmi0, FUN = length)
names(split_n) <- c("has_dm", "af_group", "pre_documented", "af_at_admission_only", "n")
write_both(split_n, "B04_preexisting_split_counts.csv", "n")

make_tvc_splitpre <- function(df, time_var, status_var) {
  d <- df
  d$g4 <- as.character(d$af_group)
  d$g4[d$af_group == "preexisting" & d$pre_documented == 1] <- "pre_documented"
  d$g4[d$af_group == "preexisting" & d$af_at_admission_only == 1] <- "af_at_admission_only"
  d$tstop0 <- d[[time_var]]; d$ev <- d[[status_var]]
  d <- d[!is.na(d$tstop0) & d$tstop0 > 0, ]
  split_noaf <- d$g4 == "NOAF" & !is.na(d$t_af_days) & d$t_af_days < d$tstop0
  a <- d; a$tstart <- 0
  a$tstop <- ifelse(split_noaf, d$t_af_days, d$tstop0)
  a$event <- ifelse(split_noaf, 0L, d$ev)
  a$exp <- ifelse(d$g4 == "NOAF", "no_AF", d$g4)
  b <- d[split_noaf, , drop = FALSE]
  if (nrow(b)) {
    b$tstart <- b$t_af_days; b$tstop <- b$tstop0; b$event <- b$ev; b$exp <- "NOAF"
    out <- rbind(a, b)
  } else out <- a
  out$exp <- factor(out$exp, levels = c("no_AF", "NOAF", "pre_documented", "af_at_admission_only"))
  out[out$tstop > out$tstart, ]
}

# ---------------- Item 3: in-hospital death competing risk + CIF figure ----------------
log("Item 3: hospital death competing risk (live discharge)")
# Cause-specific TVC Cox for death and for live discharge; Fine-Gray if survival::finegray works
# Descriptive CIF by baseline AF group (replace Fig S2)

# Cause-specific times
dmi0$etime_cr <- pmax(ifelse(is.na(dmi0$time_hosp), 1 / 24, dmi0$time_hosp), 1 / 24)
# event: 1 = hospital death, 2 = live discharge (everyone should have one)
dmi0$ev_cr <- ifelse(dmi0$status_hosp == 1, 1L, 2L)

cs_rows <- list()
# cause-specific: treat the other event as censor
for (lab in dm_labs) {
  for (k in if (isTRUE(do_mice)) 1:10 else 1) {
    dfill <- subset_dm(fill_imp(k), lab)
    dfill$etime_cr <- pmax(ifelse(is.na(dfill$time_hosp), 1 / 24, dfill$time_hosp), 1 / 24)
    dfill$ev_cr <- ifelse(dfill$status_hosp == 1, 1L, 2L)
    # death as event, discharge censored
    dfill$status_cs_death <- as.integer(dfill$ev_cr == 1)
    dfill$status_cs_dc <- as.integer(dfill$ev_cr == 2)
    # reuse make_tvc with etime
    dfill$time_cr <- dfill$etime_cr
  }
}

# Complete-case + MI pooled cause-specific and Fine-Gray
fit_cs_mi <- function(lab, evcode, model_lab) {
  coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
  for (k in 1:10) {
    dfill <- subset_dm(fill_imp(k), lab)
    dfill$time_cr <- pmax(ifelse(is.na(dfill$time_hosp), 1 / 24, dfill$time_hosp), 1 / 24)
    dfill$status_cs <- as.integer(ifelse(dfill$status_hosp == 1, 1L, 2L) == evcode)
    fit <- fit_tvc_one(dfill, "time_cr", "status_cs")
    if (is.null(fit)) next
    cf <- coef(fit); vc <- vcov(fit)
    keep <- grepl("^exp", names(cf))
    coefs[[length(coefs) + 1]] <- cf[keep]
    vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
    ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill$status_cs, na.rm = TRUE))
  }
  if (length(coefs) < 3) return(NULL)
  po <- pool_rubin(coefs, vcovs)
  po$model <- model_lab; po$subset <- lab
  po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
  po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
  po
}

cr_out <- list()
for (lab in dm_labs) {
  a <- fit_cs_mi(lab, 1, "cause_specific_hosp_death_mice")
  b <- fit_cs_mi(lab, 2, "cause_specific_live_discharge_mice")
  if (!is.null(a)) cr_out[[length(cr_out) + 1]] <- a
  if (!is.null(b)) cr_out[[length(cr_out) + 1]] <- b
}

# Fine-Gray with time-varying NOAF via survival::finegray if possible (CC, DM and all)
fg_ok <- TRUE
fg_rows <- list()
for (lab in c("DM", "all")) {
  dfill <- subset_dm(dmi0, lab)
  dfill <- dfill[stats::complete.cases(dfill[, c("time_hosp", "status_hosp", "af_group", covars)]), ]
  dfill <- drop_noaf_missing_onset(dfill, paste("FG", lab))
  dfill$etime <- pmax(dfill$time_hosp, 1 / 24)
  dfill$ev <- factor(ifelse(dfill$status_hosp == 1, "death", "discharge"),
                     levels = c("censor", "death", "discharge"))
  # counting process then finegray
  long <- make_tvc(dfill, "etime", "status_hosp")
  # rebuild event type on long: event==1 is death; if event==0 and tstop==etime and status_hosp==0, discharge
  long$ev2 <- "censor"
  long$ev2[long$event == 1] <- "death"
  # discharge occurs at etime when not death
  long$ev2[long$event == 0 & abs(long$tstop - long$etime) < 1e-8 & dfill$status_hosp[match(long$stay_id, dfill$stay_id)] == 0] <- "discharge"
  # safer: join status
  long$status_hosp <- dfill$status_hosp[match(long$stay_id, dfill$stay_id)]
  long$ev2 <- "censor"
  long$ev2[long$event == 1] <- "death"
  terminal <- long$tstop >= (dfill$etime[match(long$stay_id, dfill$stay_id)] - 1e-10)
  long$ev2[terminal & long$event == 0 & long$status_hosp == 0] <- "discharge"
  long$ev2 <- factor(long$ev2, levels = c("censor", "death", "discharge"))
  long <- restore_factors(long)
  fg_form <- as.formula(paste("Surv(tstart, tstop, ev2) ~ exp +", paste(covars, collapse = " + ")))
  fgdat <- try(survival::finegray(fg_form, data = long, id = stay_id, etype = "death"), silent = TRUE)
  if (inherits(fgdat, "try-error")) {
    log("finegray TVC failed", lab, as.character(fgdat)[1]); fg_ok <- FALSE
    next
  }
  f <- as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~ exp +", paste(covars, collapse = " + ")))
  fitfg <- try(coxph(f, data = fgdat, weights = fgwt, ties = "efron"), silent = TRUE)
  if (inherits(fitfg, "try-error")) { log("fg coxph failed", lab); next }
  ex <- extract_hr(fitfg, "finegray_tvc_cc", "hosp_death_CIF", lab,
                   n = nrow(dfill), n_event = sum(dfill$status_hosp))
  if (!is.null(ex)) {
    ex$hr_ci <- hr_fmt(ex$hr, ex$lcl, ex$ucl)
    fg_rows[[length(fg_rows) + 1]] <- ex
  }
  rm(fgdat)
}
if (length(fg_rows)) cr_out[[length(cr_out) + 1]] <- do.call(rbind, fg_rows)

# Also a simpler baseline-group Fine-Gray (no TVC) for descriptive comparison and n check
fg_base_rows <- list()
for (lab in c("DM", "all")) {
  s <- subset_dm(dmi0, lab)
  s$exp <- factor(s$af_group, levels = c("no_AF", "NOAF", "preexisting"))
  s$ftime <- pmax(s$time_hosp, 1 / 24)
  s$fstatus <- ifelse(s$status_hosp == 1, 1, 2)
  need <- c("ftime", "fstatus", "exp", covars)
  ok <- stats::complete.cases(s[, need])
  log("FG baseline complete", lab, "n_ok", sum(ok), "n_raw", nrow(s))
  if (requireNamespace("cmprsk", quietly = TRUE)) {
    X <- model.matrix(~ exp + anchor_age + sofa + sapsii + sex + charlson_ex_dm + sepsis3 +
                        cardiac_surgery + vent24 + vaso24, data = s[ok, ])[, -1, drop = FALSE]
    ok2 <- stats::complete.cases(X)
    fit <- try(cmprsk::crr(s$ftime[ok][ok2], s$fstatus[ok][ok2], X[ok2, ], failcode = 1, cencode = 0), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cf <- fit$coef; se <- sqrt(diag(fit$var)); nm <- names(cf)
      for (i in seq_along(cf)) if (grepl("^exp", nm[i])) {
        fg_base_rows[[length(fg_base_rows) + 1]] <- data.frame(
          term = sub("^exp", "", nm[i]), loghr = NA, se = se[i],
          hr = exp(cf[i]), lcl = exp(cf[i] - 1.96 * se[i]), ucl = exp(cf[i] + 1.96 * se[i]),
          p = 2 * pnorm(-abs(cf[i] / se[i])), m = NA,
          model = "finegray_crr_baseline_cc", subset = lab,
          n = sum(ok2), n_event = sum(s$fstatus[ok][ok2] == 1),
          hr_ci = hr_fmt(exp(cf[i]), exp(cf[i] - 1.96 * se[i]), exp(cf[i] + 1.96 * se[i])),
          stringsAsFactors = FALSE)
      }
    }
  }
}

item3 <- bind_fill(c(cr_out, if (length(fg_base_rows)) list(do.call(rbind, fg_base_rows))))
write_both(item3, "B03_hosp_death_competing_risk.csv", c("n", "n_event"))
log("Item 3 models done")

# CIF figure (descriptive, baseline AF group; hospital death; faceted by diabetes)
# Use cmprsk::cuminc or survival
plot_cif <- function() {
  d <- dmi0
  d$group <- factor(as.character(d$af_group), levels = c("no_AF", "NOAF", "preexisting"),
                    labels = c("No AF", "NOAF", "Pre-existing AF"))
  d$dmf <- factor(d$has_dm, levels = c(1, 0), labels = c("Diabetes", "No diabetes"))
  d$ftime <- pmax(d$time_hosp, 1 / 24)
  d$fstatus <- ifelse(d$status_hosp == 1, 1, 2)
  # truncate display at 28 days
  d$ftime28 <- pmin(d$ftime, 28)
  d$fstatus28 <- ifelse(d$ftime > 28, 0, d$fstatus)
  rows <- list()
  for (dmv in levels(d$dmf)) {
    s <- d[d$dmf == dmv, ]
    if (requireNamespace("cmprsk", quietly = TRUE)) {
      ci <- try(cmprsk::cuminc(s$ftime28, s$fstatus28, s$group, cencode = 0), silent = TRUE)
      if (inherits(ci, "try-error")) next
      # extract death curves (failcode 1): names like "No AF 1"
      nms <- names(ci)
      nms <- nms[nms != "Tests"]
      for (nm in nms) {
        if (!grepl(" 1$", nm)) next  # death
        gp <- sub(" 1$", "", nm)
        xx <- ci[[nm]]
        rows[[length(rows) + 1]] <- data.frame(dm = dmv, group = gp, time = xx$time, cif = xx$est,
                                               var = xx$var, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(NULL)
  pdat <- do.call(rbind, rows)
  pdat$group <- factor(pdat$group, levels = c("No AF", "NOAF", "Pre-existing AF"))
  # step plot
  pal <- c("No AF" = "#0072B2", "NOAF" = "#D55E00", "Pre-existing AF" = "#009E73")
  p <- ggplot(pdat, aes(x = time, y = cif, color = group)) +
    geom_step(linewidth = 0.7) +
    facet_wrap(~dm) +
    scale_color_manual(values = pal, name = NULL) +
    scale_x_continuous(limits = c(0, 28), breaks = seq(0, 28, 7)) +
    scale_y_continuous(limits = c(0, NA), labels = function(z) sprintf("%.2f", z)) +
    labs(x = "Days from ICU admission", y = "Cumulative incidence of hospital death",
         caption = "Descriptive CIF; live hospital discharge is the competing event. Follow-up displayed to day 28. Adjusted inference uses time-varying models.") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"),
          legend.position = "bottom")
  p
}
p_cif <- try(plot_cif(), silent = TRUE)
if (!inherits(p_cif, "try-error") && !is.null(p_cif)) {
  ggsave(file.path(fig_dir, "FigureS2_CIF_hospital_death.pdf"), p_cif, width = 8.5, height = 5.2, device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(fig_dir, "FigureS2_CIF_hospital_death.png"), p_cif, width = 8.5, height = 5.2, dpi = 300, device = "png")
  log("wrote CIF figure FigureS2_CIF_hospital_death")
} else log("CIF figure failed", if (inherits(p_cif, "try-error")) as.character(p_cif)[1] else "null")

# Fine-Gray n-1 note for original ICU-discharge model
n_dm_in3p <- sum(dmi0$has_dm == 1)
n_dm_fg_old <- sum(dmi0$has_dm == 1 & stats::complete.cases(dmi0[, c("icu_los_days", "died_in_icu", "af_group", "anchor_age", "sofa")]))
fg_n_note <- data.frame(
  n_dm_three_group = n_dm_in3p,
  n_dm_complete_for_original_crr_age_sofa = n_dm_fg_old,
  difference = n_dm_in3p - n_dm_fg_old,
  explanation = "Original 11_finegray_icu_discharge.csv used complete.cases on icu_los_days, died_in_icu, exp, and model.matrix(~exp+anchor_age+sofa). A difference of 1 is consistent with a single diabetic stay missing sofa (or icu_los).",
  stringsAsFactors = FALSE
)
write_both(fg_n_note, "B03_finegray_n_minus1_note.csv", c("n_dm_three_group", "n_dm_complete_for_original_crr_age_sofa", "difference"))

# ---------------- Item 4: time-zero fixes ----------------
log("Item 4: preexisting split, admission-only covars, 24h and corrected 48h landmarks, early vs late NOAF")
item4_rows <- list()

# 4a split preexisting, MI, 28-day (and hosp, 1y) DM / all
for (oc in names(outcomes)) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  for (lab in c("DM", "nonDM", "all")) {
    coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
    for (k in 1:10) {
      dfill <- subset_dm(fill_imp(k), lab)
      dfill$pre_documented <- as.integer(dfill$prior_af_icd == 1 | dfill$af_pre == 1)
      dfill$af_at_admission_only <- as.integer(as.character(dfill$af_group) == "preexisting" &
                                                 dfill$prior_af_icd == 0 & dfill$af_pre == 0)
      dfill <- drop_noaf_missing_onset(dfill, paste("splitpre", oc, lab))
      if (nrow(dfill) < 50) next
      long <- restore_factors(make_tvc_splitpre(dfill, tv, sv))
      f <- as.formula(paste("Surv(tstart, tstop, event) ~ exp +", paste(covars, collapse = " + ")))
      fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
      if (inherits(fit, "try-error")) next
      cf <- coef(fit); vc <- vcov(fit)
      keep <- grepl("^exp", names(cf))
      coefs[[length(coefs) + 1]] <- cf[keep]
      vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
      ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill[[sv]], na.rm = TRUE))
    }
    if (length(coefs) < 3) next
    po <- pool_rubin(coefs, vcovs)
    po$model <- "tvc_mice_split_preexisting"; po$outcome <- oc; po$subset <- lab
    po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
    po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
    item4_rows[[length(item4_rows) + 1]] <- po
  }
}

# 4b admission-only covariates, MI, 28-day primarily (all 3 outcomes, DM/all)
fit_mi_rhs <- function(lab, oc, rhs, model) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
  for (k in 1:10) {
    dfill <- subset_dm(fill_imp(k), lab)
    fit <- fit_tvc_one(dfill, tv, sv, rhs = rhs)
    if (is.null(fit)) next
    cf <- coef(fit); vc <- vcov(fit)
    keep <- grepl("^exp", names(cf))
    coefs[[length(coefs) + 1]] <- cf[keep]
    vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
    ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill[[sv]], na.rm = TRUE))
  }
  if (length(coefs) < 3) return(NULL)
  po <- pool_rubin(coefs, vcovs)
  po$model <- model; po$outcome <- oc; po$subset <- lab
  po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
  po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
  po
}
for (oc in names(outcomes)) {
  for (lab in c("DM", "all")) {
    a <- fit_mi_rhs(lab, oc, paste(c("exp", covars_adm), collapse = " + "), "tvc_mice_admission_covars")
    b <- fit_mi_rhs(lab, oc, paste(c("exp", covars_adm2), collapse = " + "), "tvc_mice_admission_covars_plus_csurg_sepsis")
    if (!is.null(a)) item4_rows[[length(item4_rows) + 1]] <- a
    if (!is.null(b)) item4_rows[[length(item4_rows) + 1]] <- b
  }
}

# 4c 24-h landmark: alive and still in hospital at 24h; exclude onset <=24h; NOAF if t_af>1
# 4d corrected 48-h: exclude later NOAF from the reference (do not recode them as no_AF)
landmark_fit <- function(horizon_days, mode, lab, oc) {
  # mode: "exclude_early_noaf" (24h: drop t_af<=h, late NOAF vs no_AF vs pre)
  #       "corrected" (48h: drop t_af>h from the sample; NOAF if t_af<=h)
  #       "old_style" (recode t_af>h as no_AF) -- for comparison only
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  d <- subset_dm(dmi0, lab)
  # at risk at horizon: followed past horizon and not already dead; still in hospital if hosp time
  d <- d[d[[tv]] > horizon_days, ]
  d <- d[!(d[[sv]] == 1 & d[[tv]] <= horizon_days), ]
  if (oc != "hosp_death") {
    # also require still in hospital at horizon when that time is known
    d <- d[is.na(d$time_hosp) | d$time_hosp > horizon_days, ]
  } else {
    d <- d[d$time_hosp > horizon_days, ]
  }
  if (mode == "exclude_early_noaf") {
    # 24h spec: onset <= horizon excluded; late NOAF remains NOAF (onset after landmark)
    d <- d[!(d$af_group == "NOAF" & !is.na(d$t_af_days) & d$t_af_days <= horizon_days), ]
    d$lm_group <- as.character(d$af_group)
  } else if (mode == "corrected") {
    # drop later NOAF from the reference
    d <- d[!(d$af_group == "NOAF" & (is.na(d$t_af_days) | d$t_af_days > horizon_days)), ]
    d$lm_group <- as.character(d$af_group)
  } else if (mode == "old_style") {
    d$lm_group <- as.character(d$af_group)
    d$lm_group[d$af_group == "NOAF" & (is.na(d$t_af_days) | d$t_af_days > horizon_days)] <- "no_AF"
  }
  d$lm_group <- factor(d$lm_group, levels = c("no_AF", "NOAF", "preexisting"))
  d$lm_time <- pmax(d[[tv]] - horizon_days, 1 / 24)
  d$lm_status <- d[[sv]]
  dcc <- d[stats::complete.cases(d[, c("lm_time", "lm_status", "lm_group", covars)]), ]
  if (nrow(dcc) < 50 || sum(dcc$lm_status) < 10) return(NULL)
  dcc$exp <- dcc$lm_group
  f <- as.formula(paste("Surv(lm_time, lm_status) ~ exp +", paste(covars, collapse = " + ")))
  fit <- try(coxph(f, data = dcc, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  ex <- extract_hr(fit, paste0("landmark_", horizon_days, "d_", mode, "_cc"), oc, lab,
                   n = nrow(dcc), n_event = sum(dcc$lm_status))
  if (!is.null(ex)) ex$hr_ci <- hr_fmt(ex$hr, ex$lcl, ex$ucl)
  ex
}

for (lab in dm_labs) {
  for (oc in c("death_28d", "death_1y", "hosp_death")) {
    x <- landmark_fit(1, "exclude_early_noaf", lab, oc)
    if (!is.null(x)) item4_rows[[length(item4_rows) + 1]] <- x
  }
  for (oc in c("death_28d", "death_1y")) {
    x <- landmark_fit(2, "corrected", lab, oc)
    if (!is.null(x)) item4_rows[[length(item4_rows) + 1]] <- x
    y <- landmark_fit(2, "old_style", lab, oc)
    if (!is.null(y)) item4_rows[[length(item4_rows) + 1]] <- y
  }
}

# 4e early (<=48h) vs late (>48h) NOAF as two tdc's, MI, 28-day and 1-year
make_tvc_early_late <- function(df, time_var, status_var, cut = 2) {
  d <- df
  d$tstop0 <- d[[time_var]]; d$ev <- d[[status_var]]
  d <- d[!is.na(d$tstop0) & d$tstop0 > 0, ]
  d$g <- as.character(d$af_group)
  is_noaf <- d$g == "NOAF" & !is.na(d$t_af_days) & d$t_af_days < d$tstop0
  # rows: start; optional switch at t_af
  a <- d; a$tstart <- 0
  a$tstop <- ifelse(is_noaf, d$t_af_days, d$tstop0)
  a$event <- ifelse(is_noaf, 0L, d$ev)
  a$early_noaf <- 0L; a$late_noaf <- 0L
  a$exp_pre <- as.integer(d$g == "preexisting")
  b <- d[is_noaf, , drop = FALSE]
  if (nrow(b)) {
    b$tstart <- b$t_af_days; b$tstop <- b$tstop0; b$event <- b$ev
    b$early_noaf <- as.integer(b$t_af_days <= cut)
    b$late_noaf <- as.integer(b$t_af_days > cut)
    b$exp_pre <- 0L
    out <- rbind(a, b)
  } else out <- a
  out[out$tstop > out$tstart, ]
}

for (oc in c("death_28d", "death_1y", "hosp_death")) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  for (lab in c("DM", "all")) {
    coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
    for (k in 1:10) {
      dfill <- subset_dm(fill_imp(k), lab)
      dfill <- drop_noaf_missing_onset(dfill, paste("earlylate", oc, lab))
      long <- restore_factors(make_tvc_early_late(dfill, tv, sv, cut = 2))
      f <- as.formula(paste("Surv(tstart, tstop, event) ~ early_noaf + late_noaf + exp_pre +",
                            paste(covars, collapse = " + ")))
      fit <- try(coxph(f, data = long, ties = "efron"), silent = TRUE)
      if (inherits(fit, "try-error")) next
      cf <- coef(fit); vc <- vcov(fit)
      keep <- names(cf) %in% c("early_noaf", "late_noaf", "exp_pre")
      coefs[[length(coefs) + 1]] <- cf[keep]
      vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
      ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill[[sv]], na.rm = TRUE))
    }
    if (length(coefs) < 3) next
    po <- pool_rubin(coefs, vcovs)
    po$model <- "tvc_mice_early_late48h"; po$outcome <- oc; po$subset <- lab
    po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
    po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
    item4_rows[[length(item4_rows) + 1]] <- po
  }
}

item4 <- bind_fill(item4_rows)
write_both(item4, "B04_timezero_fixes.csv", c("n", "n_event"))
log("Item 4 done, rows", nrow(item4))

# Optional strong::ok_mult_cox PHA for time-fixed corrected 48h landmark (DM, 28-day)
if (has_strong) {
  log("Item 4/7 extra: strong::ok_mult_cox PHA on corrected 48h landmark (time-fixed Cox)")
  try({
    lab <- "DM"; oc <- "death_28d"; horizon_days <- 2
    d <- subset_dm(dmi0, lab)
    d <- d[d$time_28 > horizon_days, ]
    d <- d[!(d$af_group == "NOAF" & (is.na(d$t_af_days) | d$t_af_days > horizon_days)), ]
    d$lm_group <- factor(as.character(d$af_group), levels = c("no_AF", "NOAF", "preexisting"))
    d$lm_time <- pmax(d$time_28 - horizon_days, 1 / 24)
    d$lm_status <- d$status_28
    dcc <- d[stats::complete.cases(d[, c("lm_time", "lm_status", "lm_group", covars)]), ]
    # do not print dcc
    res <- strong::ok_mult_cox(
      time = "lm_time", outcome = "lm_status", col_name_gp = "lm_group",
      covars = covars, df = dcc, num_covar = covars_num,
      op_tdc_dsp = 1, op_tdc_tt = FALSE, p_pha_global = 0.05, p_pha_covar = 0.05,
      mod_label = "landmark48h_corrected_DM_28d"
    )
    if (!is.null(res$tbl_res)) write_both(as.data.frame(res$tbl_res), "B04_strong_ok_mult_cox_lm48_tbl.csv")
    if (!is.null(res$tbl_mod)) write_both(as.data.frame(res$tbl_mod), "B04_strong_ok_mult_cox_lm48_mod.csv")
    if (!is.null(res$tbl_pha)) write_both(as.data.frame(res$tbl_pha), "B07_strong_pha_lm48.csv")
    # If tbl_pha not returned unless tt=TRUE, still try cox.zph below for TVC
    log("strong::ok_mult_cox landmark48 completed")
  }, silent = FALSE)
}

# ---------------- Item 6: ECG-available comparable subcohort ----------------
log("Item 6: ECG-available subcohort")
item6_rows <- list()
# All groups must have hospital-window ECG; ECG-negative NOAF excluded
for (lab in c("DM", "all")) {
  coefs <- list(); vcovs <- list(); ns <- c(); evs <- c(); n_excl <- c()
  for (k in 1:10) {
    dfill <- subset_dm(fill_imp(k), lab)
    dfill <- dfill[dfill$has_ecg_hosp == 1, ]
    # exclude ECG-negative NOAF
    drop <- dfill$af_group == "NOAF" & (is.na(dfill$af_hosp_ecg) | dfill$af_hosp_ecg == 0)
    n_excl <- c(n_excl, sum(drop))
    dfill <- dfill[!drop, ]
    fit <- fit_tvc_one(dfill, "time_28", "status_28")
    if (is.null(fit)) next
    cf <- coef(fit); vc <- vcov(fit)
    keep <- grepl("^exp", names(cf))
    coefs[[length(coefs) + 1]] <- cf[keep]
    vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
    ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill$status_28, na.rm = TRUE))
  }
  if (length(coefs) < 3) next
  po <- pool_rubin(coefs, vcovs)
  po$model <- "tvc_mice_ecg_available_subcohort"; po$outcome <- "death_28d"; po$subset <- lab
  po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
  po$n_noaf_ecg_neg_excluded <- round(mean(n_excl))
  po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
  item6_rows[[length(item6_rows) + 1]] <- po
}
# Descriptive: among NOAF with ECG, ECG+ vs ECG- 28-day mortality (counts)
ecg_desc <- aggregate(cbind(n = stay_id, n_death28 = status_28) ~ has_dm + I(as.character(af_group)) + has_ecg_hosp + af_hosp_ecg,
                      data = dmi0, FUN = function(z) if (identical(deparse(substitute(z)), "stay_id")) length(z) else sum(z))
# simpler explicit
ecg_desc <- {
  tmp <- dmi0
  tmp$g <- as.character(tmp$af_group)
  out <- as.data.frame(stats::xtabs(~ g + has_ecg_hosp + af_hosp_ecg + has_dm, data = tmp))
  ev <- as.data.frame(stats::xtabs(status_28 ~ g + has_ecg_hosp + af_hosp_ecg + has_dm, data = tmp))
  names(out)[names(out) == "Freq"] <- "n"
  ev$n_death28 <- ev$Freq; ev$Freq <- NULL
  merge(out, ev, by = c("g", "has_ecg_hosp", "af_hosp_ecg", "has_dm"), all = TRUE)
}
item6 <- do.call(rbind, item6_rows)
write_both(item6, "B06_ecg_available_subcohort_hr.csv", c("n", "n_event", "n_noaf_ecg_neg_excluded"))
write_both(ecg_desc, "B06_ecg_available_counts.csv", c("n", "n_death28"))
log("Item 6 done")

# ---------------- Item 7: Schoenfeld PH + strata by cardiac surgery and sepsis ----------------
log("Item 7: cox.zph and TVC strata")
zph_rows <- list()
# complete-case TVC primary models
for (oc in names(outcomes)) {
  tv <- outcomes[[oc]][1]; sv <- outcomes[[oc]][2]
  for (lab in dm_labs) {
    dfill <- subset_dm(dmi0, lab)
    dfill <- dfill[stats::complete.cases(dfill[, c(tv, sv, "af_group", covars)]), ]
    fit <- fit_tvc_one(dfill, tv, sv)
    if (is.null(fit)) next
    z <- try(cox.zph(fit, transform = "km"), silent = TRUE)
    if (inherits(z, "try-error")) { log("cox.zph failed", oc, lab); next }
    tab <- as.data.frame(z$table)
    tab$term <- rownames(z$table)
    tab$model <- "tvc_complete_case"; tab$outcome <- oc; tab$subset <- lab
    tab$n <- nrow(dfill); tab$n_event <- sum(dfill[[sv]])
    zph_rows[[length(zph_rows) + 1]] <- tab
  }
}
zph <- do.call(rbind, zph_rows)
rownames(zph) <- NULL
write_both(zph, "B07_cox_zph.csv", c("n", "n_event"))

# If exposure PH violated, time-split at 7 days (28-day DM as illustration)
split7_rows <- list()
d_dm <- subset_dm(dmi0, "DM")
d_dm <- d_dm[stats::complete.cases(d_dm[, c("time_28", "status_28", "af_group", covars)]), ]
long <- restore_factors(make_tvc(drop_noaf_missing_onset(d_dm, "split7"), "time_28", "status_28"))
# survSplit
spl <- try(survSplit(Surv(tstart, tstop, event) ~ ., data = long, cut = 7, episode = "period", id = "rowid"), silent = TRUE)
if (!inherits(spl, "try-error")) {
  spl$period <- factor(spl$period)
  f <- as.formula(paste("Surv(tstart, tstop, event) ~ exp:period + period +", paste(covars, collapse = " + ")))
  fit <- try(coxph(f, data = spl, ties = "efron"), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    sm <- summary(fit)$conf.int
    rn <- rownames(sm)
    keep <- grepl("exp", rn)
    split7_rows <- data.frame(term = rn[keep], hr = sm[keep, 1], lcl = sm[keep, 3], ucl = sm[keep, 4],
                              p = summary(fit)$coefficients[rn[keep], "Pr(>|z|)"],
                              model = "tvc_cc_survSplit_7d", outcome = "death_28d", subset = "DM",
                              n = nrow(d_dm), n_event = sum(d_dm$status_28),
                              stringsAsFactors = FALSE)
    split7_rows$hr_ci <- hr_fmt(split7_rows$hr, split7_rows$lcl, split7_rows$ucl)
    write_both(split7_rows, "B07_survsplit_7d_dm_28d.csv", c("n", "n_event"))
  }
}

# Strata: cardiac surgery and Sepsis-3, 28-day TVC, MI, DM and all
strat_rows <- list()
for (svar in c("cardiac_surgery", "sepsis3")) {
  for (lev in c(0, 1)) {
    for (lab in c("DM", "all")) {
      coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
      for (k in 1:10) {
        dfill <- subset_dm(fill_imp(k), lab)
        dfill <- dfill[dfill[[svar]] == lev, ]
        fit <- fit_tvc_one(dfill, "time_28", "status_28")
        if (is.null(fit)) next
        cf <- coef(fit); vc <- vcov(fit)
        keep <- grepl("^exp", names(cf))
        coefs[[length(coefs) + 1]] <- cf[keep]
        vcovs[[length(vcovs) + 1]] <- vc[keep, keep, drop = FALSE]
        ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill$status_28, na.rm = TRUE))
      }
      if (length(coefs) < 3) next
      po <- pool_rubin(coefs, vcovs)
      po$model <- "tvc_mice_stratum"; po$outcome <- "death_28d"; po$subset <- lab
      po$stratum_var <- svar; po$stratum_level <- lev
      po$n <- round(mean(ns)); po$n_event <- round(mean(evs))
      po$hr_ci <- hr_fmt(po$hr, po$lcl, po$ucl)
      strat_rows[[length(strat_rows) + 1]] <- po
    }
  }
}
if (length(strat_rows)) {
  write_both(do.call(rbind, strat_rows), "B07_strata_cardsurg_sepsis.csv", c("n", "n_event"))
}
log("Item 7 done")

# ---------------- Item 8: Simon-Makuch (extended KM) to replace Fig 2 ----------------
log("Item 8: Simon-Makuch curves")
sm_plot <- function(time_var, status_var, xmax, xlab, ylab, fname) {
  d <- dmi0
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
    # build step data
    # strata names
    if (is.null(sm$strata)) {
      pieces[[length(pieces) + 1]] <- data.frame(dm = dmv, group = "all", time = sm$time, surv = sm$surv,
                                                 lower = sm$lower, upper = sm$upper, stringsAsFactors = FALSE)
    } else {
      st <- as.character(sm$strata)
      st <- sub("^exp=", "", st)
      pieces[[length(pieces) + 1]] <- data.frame(dm = dmv, group = st, time = sm$time, surv = sm$surv,
                                                 lower = sm$lower, upper = sm$upper, stringsAsFactors = FALSE)
    }
  }
  pdat <- do.call(rbind, pieces)
  pdat$group <- factor(pdat$group, levels = c("No AF", "NOAF", "Pre-existing AF"))
  pal <- c("No AF" = "#0072B2", "NOAF" = "#D55E00", "Pre-existing AF" = "#009E73")
  p <- ggplot(pdat, aes(x = time, y = 1 - surv, color = group, fill = group)) +
    geom_step(linewidth = 0.7) +
    facet_wrap(~dm) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(limits = c(0, xmax), breaks = pretty(c(0, xmax), 6)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = xlab, y = ylab,
         caption = "Simon-Makuch (extended Kaplan-Meier) curves with time-varying NOAF. Follow-up starts at ICU admission; patients move from No AF to NOAF at onset. These curves are descriptive; adjusted inference uses time-varying Cox models.") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"),
          legend.position = "bottom")
  ggsave(file.path(fig_dir, paste0(fname, ".pdf")), p, width = 8.5, height = 5.2, device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(fname, ".png")), p, width = 8.5, height = 5.2, dpi = 300, device = "png")
  log("wrote", fname)
  invisible(p)
}
try(sm_plot("time_28", "status_28", 28, "Days from ICU admission", "Cumulative incidence of death (Simon-Makuch)", "Figure2_SimonMakuch_28day"), silent = FALSE)
try(sm_plot("time_1y", "status_1y", 365, "Days from ICU admission", "Cumulative incidence of death (Simon-Makuch)", "Figure2_SimonMakuch_1year"), silent = FALSE)

# ---------------- Item 5: probabilistic bias analysis, PPV ~ Beta(400, 233) ----------------
log("Item 5: QBA PPV ~ Beta(400, 233); complete-case 28-day DM; B=500")
# episensr not installed; strong has no QBA; implement as specified in the framework
qba_one <- function(dd, alpha = 400, beta = 233, scenario = c("noAF", "pre"), B = 500, seed = SEED) {
  scenario <- match.arg(scenario)
  set.seed(seed + if (scenario == "noAF") 1 else 2)
  # complete-case for speed and stability
  d0 <- dd[stats::complete.cases(dd[, c("time_28", "status_28", "af_group", covars)]), ]
  d0 <- drop_noaf_missing_onset(d0, "QBA")
  out <- rep(NA_real_, B)
  n_noaf <- sum(d0$af_group == "NOAF")
  cand0 <- which(d0$af_group == "NOAF" & !(d0$has_ecg_hosp == 1 & d0$af_hosp_ecg == 1))
  for (b in seq_len(B)) {
    ppv <- stats::rbeta(1, alpha, beta)
    n_fp <- round((1 - ppv) * n_noaf)
    n_fp <- min(n_fp, length(cand0))
    x <- d0
    if (n_fp > 0) {
      sel <- sample(cand0, n_fp)
      if (scenario == "noAF") {
        x$af_group[sel] <- "no_AF"
        x$t_af_days[sel] <- NA
      } else {
        x$af_group[sel] <- "preexisting"
        x$t_af_days[sel] <- NA
      }
      x$af_group <- factor(x$af_group, levels = c("no_AF", "NOAF", "preexisting"))
    }
    fit <- fit_tvc_one(x, "time_28", "status_28")
    if (is.null(fit)) next
    bhat <- coef(fit)["expNOAF"]
    se <- sqrt(vcov(fit)["expNOAF", "expNOAF"])
    out[b] <- bhat + stats::rnorm(1, 0, se)
  }
  ok <- out[is.finite(out)]
  data.frame(scenario = scenario, B_requested = B, B_ok = length(ok),
             median_hr = exp(stats::median(ok)),
             mean_hr = exp(mean(ok)),
             lcl = exp(stats::quantile(ok, 0.025)),
             ucl = exp(stats::quantile(ok, 0.975)),
             ppv_prior = "Beta(400, 233)",
             stringsAsFactors = FALSE)
}
ddm <- subset_dm(dmi0, "DM")
q1 <- try(qba_one(ddm, scenario = "noAF", B = 500), silent = TRUE)
q2 <- try(qba_one(ddm, scenario = "pre", B = 500), silent = TRUE)
qba <- do.call(rbind, list(
  if (!inherits(q1, "try-error")) q1 else data.frame(scenario = "noAF", error = as.character(q1)[1], stringsAsFactors = FALSE),
  if (!inherits(q2, "try-error")) q2 else data.frame(scenario = "pre", error = as.character(q2)[1], stringsAsFactors = FALSE)
))
write_both(qba, "B05_qba_ppv_beta_400_233.csv")
log("Item 5 QBA done")

# ---------------- session info ----------------
sess <- data.frame(
  item = c("date", "timezone", "R", "survival", "mice", "cmprsk", "tidycmprsk",
           "ggplot2", "duckdb", "strong", "mimicR", "mimicR430", "episensr", "seed", "SAP"),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "Pacific/Apia (UTC+13)",
    R.version.string,
    as.character(packageVersion("survival")),
    as.character(packageVersion("mice")),
    if (requireNamespace("cmprsk", quietly = TRUE)) as.character(packageVersion("cmprsk")) else "not installed",
    if (requireNamespace("tidycmprsk", quietly = TRUE)) as.character(packageVersion("tidycmprsk")) else "not installed",
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("duckdb")),
    if (has_strong) as.character(packageVersion("strong")) else "not loaded",
    if (requireNamespace("mimicR", quietly = TRUE)) as.character(packageVersion("mimicR")) else "not installed",
    if (requireNamespace("mimicR430", quietly = TRUE)) as.character(packageVersion("mimicR430")) else "not installed",
    if (requireNamespace("episensr", quietly = TRUE)) as.character(packageVersion("episensr")) else "not installed (QBA implemented in base/survival)",
    as.character(SEED),
    "docs/SAP_v1.md; add-on does not modify locked results_final"
  ),
  stringsAsFactors = FALSE
)
write_both(sess, "B00_session_info.csv")

# disconnect duckdb if still open
try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
log("DONE B_addon")
