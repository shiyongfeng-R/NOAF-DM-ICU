#!/usr/bin/env Rscript
# SAP v1.0 main analysis. Aggregates only written to disk.
# Usage (from the repository root, after setting environment variables; see README.md):
#   Rscript --vanilla scripts/01_run_main.R
#   nohup Rscript --vanilla scripts/01_run_main.R > logs/01_run_main.out 2>&1 &
source(if (file.exists("scripts/00_helpers.R")) "scripts/00_helpers.R" else "00_helpers.R")
set.seed(SEED)
log("START SAP v1.0 main analysis; seed", SEED)

catalog <- file.path(root, "catalog/coredata.duckdb")
mim_path <- mim_duckdb
if (!nzchar(root) || !file.exists(catalog)) stop("NOAF_ROOT must point to a warehouse root that contains catalog/coredata.duckdb (see README.md).")
if (!nzchar(mim_path) || !file.exists(mim_path)) stop("Set NOAF_MIMIC_DUCKDB to an existing mimic.duckdb file (READ_ONLY; see README.md).")

con <- dbConnect(duckdb::duckdb(), dbdir=":memory:")
on.exit(try(dbDisconnect(con, shutdown=TRUE), silent=TRUE), add=TRUE)
dbExecute(con, sprintf("ATTACH '%s' AS wh (READ_ONLY)", gsub("'", "''", catalog, fixed=TRUE)))
dbExecute(con, sprintf("ATTACH '%s' AS mimic (READ_ONLY)", gsub("'", "''", mim_path, fixed=TRUE)))
log("attached catalog + mimic.duckdb READ_ONLY")

# Window functions use a deterministic ORDER BY (time, identifiers, value) so tied timestamps
# cannot change first-rhythm / first-event assignment between runs (SAP v1.0 implementation fix).
# ---------------- ECG machine AF flags (same regex as feasibility validation) ----------------
log("ECG flags")
inv_mm <- dbGetQuery(con, "DESCRIBE SELECT * FROM wh.mimic_ecg.machine_measurements")
rep_cols <- inv_mm$column_name[grepl("^report(_|$)", inv_mm$column_name, ignore.case=TRUE)]
concat_sql <- paste(sprintf("CAST(%s AS VARCHAR)", paste0('"', gsub('"','""',rep_cols), '"')), collapse=", ")
dbExecute(con, sprintf("
CREATE OR REPLACE TEMP TABLE ecg_raw AS
SELECT subject_id, study_id, CAST(ecg_time AS TIMESTAMP) AS ecg_time,
       lower(trim(concat_ws(' ', %s))) AS report_text
FROM wh.mimic_ecg.machine_measurements
", concat_sql))

clean_sql <- paste0(
  "regexp_replace(regexp_replace(regexp_replace(regexp_replace(regexp_replace(regexp_replace(regexp_replace(",
  "lower(coalesce(report_text,'')), ",
  "'no longer([[:space:]]+in){0,1}[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g'), ",
  "'no[[:space:]]+(longer[[:space:]]+){0,1}evidence[[:space:]]+(of|for)[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g'), ",
  "'ruled[[:space:]]+out[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g'), ",
  "'(absence[[:space:]]+of|negative[[:space:]]+for|without|denies|not)[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g'), ",
  "'no[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g'), ",
  "'(history|hx|past[[:space:]]+history)[[:space:]]+of[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g'), ",
  "'(resolved|has[[:space:]]+resolved)[[:space:]]+(atrial[[:space:]]+fibrillation|atrial[[:space:]]+flutter|a-{0,1}fib(rillation){0,1}|a-{0,1}flutter)', ' ', 'g')"
)
af_re <- "atrial[[:space:]]+fibrillat|a-fib|[[:space:]]afib[[:space:]]|^afib[[:space:]]|[[:space:]]afib$|^afib$|atrial[[:space:]]+fib([[:space:]]|$)"
afl_re <- "atrial[[:space:]]+flutter|a-flutter|aflutter"
poss_re <- paste0(
  "(possible|consider|probable|questionable|cannot[[:space:]]+rule[[:space:]]+out).{0,50}",
  "(atrial[[:space:]]+fibrillat|a-fib|afib|atrial[[:space:]]+fib|atrial[[:space:]]+flutter|a-flutter|aflutter)",
  "|(atrial[[:space:]]+fibrillat|a-fib|afib|atrial[[:space:]]+fib|atrial[[:space:]]+flutter|a-flutter|aflutter)",
  ".{0,40}(possible|cannot[[:space:]]+be[[:space:]]+ruled[[:space:]]+out)"
)
dbExecute(con, sprintf("
CREATE OR REPLACE TEMP TABLE ecg AS
SELECT subject_id, study_id, ecg_time,
  CASE WHEN regexp_matches(%s, '%s') OR regexp_matches(%s, '%s') THEN 1 ELSE 0 END AS has_af_any,
  CASE WHEN (regexp_matches(%s, '%s') OR regexp_matches(%s, '%s'))
            AND regexp_matches(%s, '%s')=FALSE THEN 1 ELSE 0 END AS has_af_definite
FROM ecg_raw
", clean_sql, af_re, clean_sql, afl_re,
   clean_sql, af_re, clean_sql, afl_re, clean_sql, poss_re))
try(dbExecute(con, "DROP TABLE ecg_raw"), silent=TRUE)
log("ECG flagged")

# ---------------- First ICU adult + ICD ----------------
log("first ICU adult")
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE first_icu AS
WITH first_icu AS (
  SELECT
    i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.los,
    p.anchor_age, p.anchor_year_group, p.gender, p.dod,
    a.admittime, a.dischtime, a.deathtime, a.hospital_expire_flag,
    a.admission_type, a.race, a.insurance
  FROM wh.mimic_iv.icustays i
  INNER JOIN wh.mimic_iv.patients p USING (subject_id)
  INNER JOIN wh.mimic_iv.admissions a USING (subject_id, hadm_id)
  WHERE p.anchor_age >= 18
  QUALIFY ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.intime, i.stay_id, i.hadm_id) = 1
),
dx AS (
  SELECT hadm_id,
    MAX(CASE WHEN (icd_version=9 AND replace(icd_code,'.','')='42731')
              OR (icd_version=10 AND regexp_matches(upper(replace(trim(icd_code),'.','')), '^I48'))
            THEN 1 ELSE 0 END) AS index_af_icd,
    MAX(CASE WHEN (icd_version=9 AND replace(icd_code,'.','') LIKE '250%')
              OR (icd_version=10 AND regexp_matches(upper(replace(trim(icd_code),'.','')), '^E1[0-4]'))
            THEN 1 ELSE 0 END) AS has_dm_icd,
    MAX(CASE WHEN (icd_version=9 AND regexp_matches(replace(icd_code,'.',''), '^428'))
              OR (icd_version=10 AND regexp_matches(upper(replace(trim(icd_code),'.','')), '^I50|^I110|^I130|^I132|^I255'))
            THEN 1 ELSE 0 END) AS has_hf,
    MAX(CASE WHEN (icd_version=9 AND regexp_matches(replace(icd_code,'.',''), '^40[1-5]'))
              OR (icd_version=10 AND regexp_matches(upper(replace(trim(icd_code),'.','')), '^I1[0-5]'))
            THEN 1 ELSE 0 END) AS has_htn,
    MAX(CASE WHEN (icd_version=9 AND regexp_matches(replace(icd_code,'.',''), '^585|^403|^404|^V451|^V56'))
              OR (icd_version=10 AND regexp_matches(upper(replace(trim(icd_code),'.','')), '^N18|^N19|^I120|^I131|^Z992|^Z49'))
            THEN 1 ELSE 0 END) AS has_ckd,
    MAX(CASE WHEN (icd_version=9 AND regexp_matches(replace(icd_code,'.',''), '^43[0-4]|^436'))
              OR (icd_version=10 AND regexp_matches(upper(replace(trim(icd_code),'.','')), '^I6[0-4]'))
            THEN 1 ELSE 0 END) AS has_stroke
  FROM wh.mimic_iv.diagnoses_icd
  GROUP BY hadm_id
)
SELECT
  f.*,
  coalesce(d.index_af_icd,0) AS index_af_icd,
  coalesce(d.has_dm_icd,0) AS has_dm_icd,
  coalesce(d.has_hf,0) AS has_hf,
  coalesce(d.has_htn,0) AS has_htn,
  coalesce(d.has_ckd,0) AS has_ckd,
  coalesce(d.has_stroke,0) AS has_stroke,
  CASE WHEN f.anchor_year_group IN ('2008 - 2010','2011 - 2013','2014 - 2016','2017 - 2019') THEN 1 ELSE 0 END AS year_primary,
  CASE WHEN f.hospital_expire_flag=1 THEN 1 ELSE 0 END AS hosp_death,
  CASE
    WHEN f.dod IS NOT NULL AND f.dod <= CAST(f.intime AS DATE) + 28 THEN 1
    WHEN f.hospital_expire_flag=1 AND CAST(coalesce(f.deathtime, f.dischtime) AS DATE) <= CAST(f.intime AS DATE) + 28 THEN 1
    ELSE 0
  END AS death_28d,
  CASE
    WHEN f.dod IS NOT NULL AND f.dod <= CAST(f.intime AS DATE) + 365 THEN 1
    WHEN f.hospital_expire_flag=1 AND CAST(coalesce(f.deathtime, f.dischtime) AS DATE) <= CAST(f.intime AS DATE) + 365 THEN 1
    ELSE 0
  END AS death_1y,
  CASE WHEN f.dod IS NOT NULL THEN 1 ELSE 0 END AS has_dod,
  CASE
    WHEN f.dod IS NOT NULL THEN date_diff('day', CAST(f.intime AS DATE), f.dod)
    WHEN f.hospital_expire_flag=1 THEN date_diff('day', CAST(f.intime AS DATE), CAST(coalesce(f.deathtime, f.dischtime) AS DATE))
    ELSE NULL
  END AS days_to_death,
  date_diff('day', CAST(f.intime AS DATE), CAST(f.dischtime AS DATE)) AS days_to_discharge,
  date_diff('hour', f.intime, f.outtime)/24.0 AS icu_los_days
FROM first_icu f
LEFT JOIN dx d USING (hadm_id)
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE prior_af AS
SELECT DISTINCT i.subject_id
FROM first_icu i
INNER JOIN wh.mimic_iv.diagnoses_icd d ON d.subject_id=i.subject_id
INNER JOIN wh.mimic_iv.admissions a ON a.hadm_id=d.hadm_id AND a.subject_id=d.subject_id
WHERE a.admittime < i.admittime
  AND (
    (d.icd_version=9 AND replace(d.icd_code,'.','')='42731')
    OR (d.icd_version=10 AND regexp_matches(upper(replace(trim(d.icd_code),'.','')), '^I48'))
  )
")

# ---------------- Rhythm ----------------
log("Heart Rhythm 220048")
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE rhythm AS
SELECT
  c.stay_id,
  e.charttime, e.storetime, e.value, e.valuenum, e.caregiver_id,
  CASE WHEN regexp_matches(lower(coalesce(e.value,'')),
            'atrial fib|a-fib|afib|(^|[^a-z])af([^a-z]|$)|a flut|flutter') THEN 1 ELSE 0 END AS is_af,
  CASE WHEN regexp_matches(lower(coalesce(e.value,'')), 'sinus') THEN 1 ELSE 0 END AS is_sinus
FROM first_icu c
INNER JOIN wh.mimic_iv.chartevents e ON e.stay_id=c.stay_id
WHERE e.itemid=220048
  AND e.charttime >= c.intime - INTERVAL 1 HOUR
  AND e.charttime <= c.outtime
")
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE rhythm_stay AS
WITH seq AS (
  SELECT stay_id, charttime, storetime, value, valuenum, caregiver_id, is_af, is_sinus,
         ROW_NUMBER() OVER (PARTITION BY stay_id ORDER BY charttime, storetime, value, valuenum, caregiver_id) AS rn,
         LAG(is_af) OVER (PARTITION BY stay_id ORDER BY charttime, storetime, value, valuenum, caregiver_id) AS prev_af
  FROM rhythm
),
runs AS (
  SELECT stay_id, charttime, is_af,
         SUM(CASE WHEN is_af=1 AND coalesce(prev_af,0)=0 THEN 1 ELSE 0 END)
           OVER (PARTITION BY stay_id ORDER BY charttime, storetime, value, valuenum, caregiver_id) AS run_id
  FROM seq
),
run_len AS (
  SELECT stay_id, run_id,
         COUNT(*) AS n_in_run,
         date_diff('second', MIN(charttime), MAX(charttime))/3600.0 AS run_hours
  FROM runs
  WHERE is_af=1 AND run_id>0
  GROUP BY stay_id, run_id
)
SELECT
  s.stay_id,
  MAX(CASE WHEN s.rn=1 THEN s.value END) AS first_value,
  MIN(s.charttime) FILTER (WHERE s.rn=1) AS first_rhythm_time,
  MAX(s.is_af) FILTER (WHERE s.rn=1) AS first_rhythm_is_af,
  MAX(s.is_sinus) FILTER (WHERE s.rn=1) AS first_is_sinus,
  MIN(s.charttime) FILTER (WHERE s.is_af=1) AS first_af_time,
  COUNT(*) AS n_rhythm,
  COUNT(*) FILTER (WHERE s.is_af=1) AS n_af,
  coalesce(MAX(r.n_in_run),0) AS max_af_run_n,
  coalesce(MAX(r.run_hours),0) AS max_af_run_hours
FROM seq s
LEFT JOIN run_len r USING (stay_id)
GROUP BY s.stay_id
")

# ---------------- Exposure (SAP v1) ----------------
log("exposure classification SAP v1")
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE pre_ecg AS
SELECT c.stay_id,
  MAX(CASE WHEN e.ecg_time < c.intime AND e.has_af_any=1 THEN 1 ELSE 0 END) AS af_pre,
  MAX(CASE WHEN e.ecg_time < c.intime THEN 1 ELSE 0 END) AS has_ecg_pre,
  MAX(CASE WHEN e.ecg_time >= c.intime AND e.ecg_time <= c.outtime AND e.has_af_any=1 THEN 1 ELSE 0 END) AS af_icu_ecg,
  MAX(CASE WHEN e.ecg_time >= c.intime AND e.ecg_time <= c.outtime THEN 1 ELSE 0 END) AS has_ecg_icu,
  MAX(CASE WHEN e.ecg_time >= c.intime AND e.ecg_time <= c.dischtime AND e.has_af_any=1 THEN 1 ELSE 0 END) AS af_hosp_ecg,
  MAX(CASE WHEN e.ecg_time >= c.intime AND e.ecg_time <= c.dischtime THEN 1 ELSE 0 END) AS has_ecg_hosp
FROM first_icu c
LEFT JOIN ecg e ON e.subject_id=c.subject_id
GROUP BY c.stay_id
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE exposed AS
SELECT
  c.*,
  CASE WHEN p.subject_id IS NOT NULL THEN 1 ELSE 0 END AS prior_af_icd,
  CASE WHEN r.stay_id IS NOT NULL THEN 1 ELSE 0 END AS has_rhythm,
  coalesce(r.first_rhythm_is_af,0) AS first_rhythm_is_af,
  coalesce(r.first_is_sinus,0) AS first_is_sinus,
  r.first_af_time, r.first_rhythm_time, r.first_value,
  coalesce(r.n_af,0) AS n_af,
  coalesce(r.max_af_run_n,0) AS max_af_run_n,
  coalesce(r.max_af_run_hours,0) AS max_af_run_hours,
  CASE WHEN r.first_af_time IS NOT NULL AND r.first_af_time <= c.intime + INTERVAL 2 HOUR THEN 1 ELSE 0 END AS early_af_2h,
  coalesce(e.af_pre,0) AS af_pre,
  coalesce(e.has_ecg_pre,0) AS has_ecg_pre,
  coalesce(e.af_icu_ecg,0) AS af_icu_ecg,
  coalesce(e.has_ecg_icu,0) AS has_ecg_icu,
  coalesce(e.af_hosp_ecg,0) AS af_hosp_ecg,
  coalesce(e.has_ecg_hosp,0) AS has_ecg_hosp,
  CASE WHEN r.first_af_time IS NOT NULL
       THEN date_diff('second', c.intime, r.first_af_time)/86400.0
       ELSE NULL END AS t_af_days
FROM first_icu c
LEFT JOIN prior_af p USING (subject_id)
LEFT JOIN rhythm_stay r USING (stay_id)
LEFT JOIN pre_ecg e USING (stay_id)
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE classed AS
SELECT
  *,
  CASE WHEN prior_af_icd=1 OR af_pre=1 OR early_af_2h=1 THEN 1 ELSE 0 END AS is_preexisting,
  CASE WHEN prior_af_icd=0 AND af_pre=0 AND early_af_2h=0
            AND has_rhythm=1 AND first_is_sinus=1 AND first_rhythm_is_af=0 AND n_af>0
       THEN 1 ELSE 0 END AS is_noaf,
  CASE WHEN prior_af_icd=0 AND af_pre=0 AND early_af_2h=0
            AND has_rhythm=1 AND n_af=0 AND index_af_icd=0
       THEN 1 ELSE 0 END AS is_no_af,
  CASE WHEN index_af_icd=1 AND prior_af_icd=0 AND af_pre=0 AND early_af_2h=0
            AND has_rhythm=1 AND n_af=0
       THEN 1 ELSE 0 END AS is_grayzone,
  -- sensitivity definitions
  CASE WHEN prior_af_icd=0 AND af_pre=0 AND early_af_2h=0
            AND has_rhythm=1 AND first_rhythm_is_af=0 AND n_af>0
       THEN 1 ELSE 0 END AS is_noaf_broad,
  CASE WHEN prior_af_icd=0 AND af_pre=0 AND early_af_2h=0
            AND has_rhythm=1 AND first_is_sinus=1 AND first_rhythm_is_af=0 AND n_af>0
            AND max_af_run_n>=2
       THEN 1 ELSE 0 END AS is_noaf_2consec,
  CASE WHEN prior_af_icd=0 AND af_pre=0 AND early_af_2h=0
            AND has_rhythm=1 AND first_is_sinus=1 AND first_rhythm_is_af=0 AND n_af>0
            AND max_af_run_hours>=1
       THEN 1 ELSE 0 END AS is_noaf_1h,
  CASE WHEN prior_af_icd=1 OR early_af_2h=1 THEN 1 ELSE 0 END AS is_preexisting_noecg
FROM exposed
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE classed2 AS
SELECT *,
  CASE
    WHEN is_preexisting=1 THEN 'preexisting'
    WHEN is_noaf=1 THEN 'NOAF'
    WHEN is_no_af=1 THEN 'no_AF'
    WHEN is_grayzone=1 THEN 'grayzone'
    ELSE 'indeterminate'
  END AS af_group,
  CASE
    WHEN is_preexisting_noecg=1 THEN 'preexisting'
    WHEN prior_af_icd=0 AND early_af_2h=0 AND has_rhythm=1 AND first_is_sinus=1 AND first_rhythm_is_af=0 AND n_af>0 THEN 'NOAF'
    WHEN prior_af_icd=0 AND early_af_2h=0 AND has_rhythm=1 AND n_af=0 AND index_af_icd=0 THEN 'no_AF'
    WHEN index_af_icd=1 AND prior_af_icd=0 AND early_af_2h=0 AND has_rhythm=1 AND n_af=0 THEN 'grayzone'
    ELSE 'indeterminate'
  END AS af_group_noecg
FROM classed
")
log("classified", dbGetQuery(con,"SELECT COUNT(*) n FROM classed2")$n)

if (identical(Sys.getenv("NOAF_STOP_AFTER_CLASS"), "1")) {
  stab <- Sys.getenv("NOAF_STAB_DIR")
  if (!nzchar(stab)) stab <- file.path(proj, "logs", "stab_tmp")
  dir.create(stab, recursive=TRUE, showWarnings=FALSE)
  flow <- dbGetQuery(con, "
SELECT * FROM (
  SELECT '01_adult_first_icu' AS step, COUNT(*) AS n, COUNT(*) FILTER (WHERE has_dm_icd=1) AS n_dm FROM first_icu
  UNION ALL
  SELECT '02_year_2008_2019', COUNT(*) FILTER (WHERE year_primary=1), COUNT(*) FILTER (WHERE year_primary=1 AND has_dm_icd=1) FROM first_icu
  UNION ALL
  SELECT '03_year_2020_2022', COUNT(*) FILTER (WHERE year_primary=0), COUNT(*) FILTER (WHERE year_primary=0 AND has_dm_icd=1) FROM first_icu
  UNION ALL
  SELECT '04_has_rhythm_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND has_rhythm=1), COUNT(*) FILTER (WHERE year_primary=1 AND has_rhythm=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '05_preexisting_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '06_NOAF_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_noaf=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_noaf=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '07_no_AF_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_no_af=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_no_af=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '08_grayzone_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_grayzone=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_grayzone=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '09_indeterminate_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting=0 AND is_noaf=0 AND is_no_af=0 AND is_grayzone=0),
         COUNT(*) FILTER (WHERE year_primary=1 AND has_dm_icd=1 AND is_preexisting=0 AND is_noaf=0 AND is_no_af=0 AND is_grayzone=0) FROM classed2
  UNION ALL
  SELECT '10_three_group_primary_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting+is_noaf+is_no_af=1),
         COUNT(*) FILTER (WHERE year_primary=1 AND has_dm_icd=1 AND is_preexisting+is_noaf+is_no_af=1) FROM classed2
)
ORDER BY step
")
  exp_n <- dbGetQuery(con, "
SELECT year_primary, has_dm_icd AS has_dm, af_group, COUNT(*) AS n,
       COUNT(*) FILTER (WHERE hosp_death=1) AS n_hosp_death,
       COUNT(*) FILTER (WHERE death_28d=1) AS n_death_28d,
       COUNT(*) FILTER (WHERE death_1y=1) AS n_death_1y
FROM classed2
GROUP BY 1,2,3
ORDER BY 1 DESC, 2 DESC, 3
")
  write.csv(flow, file.path(stab, "01_flow_counts.csv"), row.names=FALSE)
  write.csv(exp_n, file.path(stab, "02_group_counts_mortality.csv"), row.names=FALSE)
  log("STOP after classification; wrote", stab)
  try(dbDisconnect(con, shutdown=TRUE), silent=TRUE)
  quit(save="no", status=0)
}

if (identical(Sys.getenv("NOAF_PUB_KM"), "1")) {
  log("PUB KM from classed2 (no individual files written)")
  dat <- dbGetQuery(con, "
    SELECT stay_id, has_dm_icd AS has_dm, af_group, year_primary,
           hosp_death, death_28d, death_1y, days_to_death, days_to_discharge, icu_los_days
    FROM classed2
    WHERE year_primary=1 AND af_group IN ('no_AF','NOAF','preexisting')
  ")
  safe_time <- function(x, default=28) pmax(ifelse(is.na(x), default, x), 1/24)
  dat$time_hosp <- safe_time(dat$days_to_discharge, 1)
  dat$status_hosp <- as.integer(dat$hosp_death==1)
  dat$time_28 <- pmin(safe_time(ifelse(dat$death_28d==1, dat$days_to_death, 28), 28), 28)
  dat$status_28 <- as.integer(dat$death_28d==1)
  dat$time_1y <- pmin(safe_time(ifelse(dat$death_1y==1, dat$days_to_death, 365), 365), 365)
  dat$status_1y <- as.integer(dat$death_1y==1)
  dat$af_group <- factor(dat$af_group, levels=c("no_AF","NOAF","preexisting"),
                         labels=c("No AF","NOAF","Pre-existing AF"))
  dat$dm_lab <- factor(dat$has_dm, levels=c(0,1), labels=c("Without diabetes","With diabetes"))
  pub_fig <- file.path(proj, "results_final", "figures")
  dir.create(pub_fig, recursive=TRUE, showWarnings=FALSE)
  save_fig_pdf <- function(p, name, width=7.4, height=5.2) {
    pdfp <- file.path(pub_fig, paste0(name, ".pdf"))
    pngp <- file.path(pub_fig, paste0(name, ".png"))
    if (file.exists(pdfp)) unlink(pdfp)
    if (file.exists(pngp)) unlink(pngp)
    ggsave(pdfp, p, width=width, height=height, device=grDevices::pdf, useDingbats=FALSE)
    ggsave(pngp, p, width=width, height=height, dpi=300, device="png")
    log("pub figure", name)
  }
  pal <- c("No AF"="#0072B2", "NOAF"="#D55E00", "Pre-existing AF"="#009E73")
  # manual ggplot KM without survminer
  km_df <- function(timevar, statusvar) {
    out <- list()
    for (dm in levels(dat$dm_lab)) for (g in levels(dat$af_group)) {
      s <- dat[dat$dm_lab==dm & dat$af_group==g, ]
      fit <- survfit(Surv(s[[timevar]], s[[statusvar]]) ~ 1)
      sm <- summary(fit)
      if (!length(sm$time)) next
      dd <- data.frame(time=c(0, sm$time), surv=c(1, sm$surv),
                       lcl=c(1, sm$lower), ucl=c(1, sm$upper),
                       af_group=g, dm_lab=dm, stringsAsFactors=FALSE)
      out[[length(out)+1]] <- dd
    }
    do.call(rbind, out)
  }
  theme_km <- theme_bw(base_size=10, base_family="sans") +
    theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
          strip.background=element_rect(fill="white", colour="black"),
          legend.position="bottom", legend.title=element_blank())
  make_km <- function(df, xlab, title, xmax) {
    df <- df[df$time<=xmax, ]
    ggplot(df, aes(time, surv, colour=af_group, fill=af_group)) +
      geom_ribbon(aes(ymin=pmax(0,lcl), ymax=pmin(1,ucl)), alpha=0.12, colour=NA) +
      geom_step(linewidth=0.55) +
      facet_wrap(~dm_lab, ncol=2) +
      scale_colour_manual(values=pal) + scale_fill_manual(values=pal) +
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2), expand=c(0.01,0.01)) +
      scale_x_continuous(limits=c(0,xmax), expand=c(0.01,0.01)) +
      labs(x=xlab, y="Survival probability", title=title,
           caption="Descriptive Kaplan-Meier curves. Groups are baseline labels. Primary inference uses time-varying Cox models.") +
      theme_km
  }
  save_fig_pdf(make_km(km_df("time_28","status_28"), "Days from ICU admission", "28-day all-cause survival", 28), "Figure2_KM_28day", 8.2, 4.6)
  save_fig_pdf(make_km(km_df("time_1y","status_1y"), "Days from ICU admission", "One-year all-cause survival", 365), "Figure2_KM_1year", 8.2, 4.6)
  save_fig_pdf(make_km(km_df("time_hosp","status_hosp"), "Days from ICU admission", "In-hospital survival (censored at discharge)", 28), "FigureS_KM_hospital_28day", 8.2, 4.6)
  try(dbDisconnect(con, shutdown=TRUE), silent=TRUE)
  log("PUB KM done")
  quit(save="no", status=0)
}

# ---------------- Covariates ----------------
log("covariates: cardiac surgery, charlson, severity, vent, vaso, labs")
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE cardiac AS
SELECT c.stay_id,
  MAX(CASE WHEN s.curr_service='CSURG' THEN 1 ELSE 0 END) AS csurg_service,
  MAX(CASE WHEN pr.icd_version=9 AND regexp_matches(replace(pr.icd_code,'.',''), '^361|^35[1-9]') THEN 1 ELSE 0 END) AS cardiac_proc9,
  MAX(CASE WHEN pr.icd_version=10 AND regexp_matches(upper(replace(trim(pr.icd_code),'.','')), '^021|^02[RQS]') THEN 1 ELSE 0 END) AS cardiac_proc10
FROM first_icu c
LEFT JOIN wh.mimic_iv.services s ON s.hadm_id=c.hadm_id
LEFT JOIN wh.mimic_iv.procedures_icd pr ON pr.hadm_id=c.hadm_id
GROUP BY c.stay_id
")

# Charlson (Quan): prefixes on codes without dots. Diabetes component stored separately.
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE charlson AS
WITH d AS (
  SELECT c.stay_id, c.hadm_id,
    upper(replace(trim(dx.icd_code),'.','')) AS code,
    dx.icd_version
  FROM first_icu c
  INNER JOIN wh.mimic_iv.diagnoses_icd dx ON dx.hadm_id=c.hadm_id
),
flag AS (
  SELECT stay_id,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^I21|^I22|^I252'))
              OR (icd_version=9 AND regexp_matches(code,'^410|^412')) THEN 1 ELSE 0 END) AS mi,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^I099|^I110|^I130|^I132|^I255|^I420|^I425|^I426|^I427|^I428|^I429|^I43|^I50'))
              OR (icd_version=9 AND regexp_matches(code,'^39891|^40201|^40211|^40291|^40401|^40403|^40411|^40413|^40491|^40493|^4254|^4255|^4257|^4258|^4259|^428')) THEN 1 ELSE 0 END) AS chf,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^I70|^I71|^I731|^I738|^I739|^I771|^I790|^I792|^K551|^K558|^K559|^Z958|^Z959'))
              OR (icd_version=9 AND regexp_matches(code,'^0930|^4373|^440|^441|^4431|^4432|^4438|^4439|^4471|^5571|^5579|^V434')) THEN 1 ELSE 0 END) AS pvd,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^G45|^G46|^H340|^I6[0-9]'))
              OR (icd_version=9 AND regexp_matches(code,'^36234|^43[0-8]')) THEN 1 ELSE 0 END) AS cevd,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^F0[0-3]|^F051|^G30|^G311'))
              OR (icd_version=9 AND regexp_matches(code,'^290|^2941|^3312')) THEN 1 ELSE 0 END) AS dem,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^I278|^I279|^J4[0-7]|^J6[0-7]|^J684|^J701|^J703'))
              OR (icd_version=9 AND regexp_matches(code,'^4168|^4169|^490|^491|^492|^493|^494|^495|^496|^497|^498|^499|^500|^501|^502|^503|^504|^505|^5064|^5081|^5088')) THEN 1 ELSE 0 END) AS cpd,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^M05|^M06|^M315|^M32|^M33|^M34|^M351|^M353|^M360'))
              OR (icd_version=9 AND regexp_matches(code,'^4465|^7100|^7101|^7102|^7103|^7104|^7140|^7141|^7142|^7148|^725')) THEN 1 ELSE 0 END) AS rheum,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^K2[5-8]'))
              OR (icd_version=9 AND regexp_matches(code,'^531|^532|^533|^534')) THEN 1 ELSE 0 END) AS pud,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^B18|^K700|^K701|^K702|^K703|^K709|^K713|^K714|^K715|^K717|^K73|^K74|^K760|^K762|^K763|^K764|^K768|^K769|^Z944'))
              OR (icd_version=9 AND regexp_matches(code,'^07022|^07023|^07032|^07033|^07044|^07054|^0706|^0709|^570|^571|^5733|^5734|^5738|^5739|^V427')) THEN 1 ELSE 0 END) AS mld,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^E1[0-4]'))
              OR (icd_version=9 AND regexp_matches(code,'^2500|^2501|^2502|^2503|^2507|^2508|^2509')) THEN 1 ELSE 0 END) AS diab,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^E1[0-4][2-57]'))
              OR (icd_version=9 AND regexp_matches(code,'^2504|^2505|^2506')) THEN 1 ELSE 0 END) AS diab_c,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^G041|^G114|^G801|^G802|^G81|^G82|^G830|^G831|^G832|^G833|^G834|^G839'))
              OR (icd_version=9 AND regexp_matches(code,'^3341|^342|^343|^3440|^3441|^3442|^3443|^3444|^3445|^3446|^3449')) THEN 1 ELSE 0 END) AS hp,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^I120|^I131|^N032|^N033|^N034|^N035|^N036|^N037|^N052|^N053|^N054|^N055|^N056|^N057|^N18|^N19|^N250|^Z49|^Z940|^Z992'))
              OR (icd_version=9 AND regexp_matches(code,'^40301|^40311|^40391|^40402|^40403|^40412|^40413|^40492|^40493|^582|^5830|^5831|^5832|^5834|^5836|^5837|^585|^586|^5880|^V420|^V451|^V56')) THEN 1 ELSE 0 END) AS renal,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^C0|^C1|^C2[0-6]|^C3[0-4]|^C37|^C38|^C39|^C4[0-1]|^C43|^C4[5-9]|^C5|^C6|^C7[0-6]|^C81|^C82|^C83|^C84|^C85|^C88|^C9[0-7]'))
              OR (icd_version=9 AND regexp_matches(code,'^14|^15|^16|^17|^18|^19|^200|^201|^202|^203|^204|^205|^206|^207|^208')) THEN 1 ELSE 0 END) AS canc,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^I850|^I859|^I864|^I982|^K704|^K711|^K721|^K729|^K765|^K766|^K767'))
              OR (icd_version=9 AND regexp_matches(code,'^4560|^4561|^4562|^5722|^5723|^5724|^5728')) THEN 1 ELSE 0 END) AS msld,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^C7[7-9]|^C80'))
              OR (icd_version=9 AND regexp_matches(code,'^19[6-9]')) THEN 1 ELSE 0 END) AS met,
    MAX(CASE WHEN (icd_version=10 AND regexp_matches(code,'^B2[0-4]'))
              OR (icd_version=9 AND regexp_matches(code,'^042|^043|^044')) THEN 1 ELSE 0 END) AS aids
  FROM d
  GROUP BY stay_id
)
SELECT stay_id,
  mi+chf+pvd+cevd+dem+cpd+rheum+pud
  + CASE WHEN msld=1 THEN 0 ELSE mld END
  + CASE WHEN diab_c=1 THEN 0 ELSE diab END
  + 2*diab_c + 2*hp + 2*renal
  + CASE WHEN met=1 THEN 0 ELSE 2*canc END
  + 3*msld + 6*met + 6*aids AS charlson,
  mi+chf+pvd+cevd+dem+cpd+rheum+pud
  + CASE WHEN msld=1 THEN 0 ELSE mld END
  + 2*hp + 2*renal
  + CASE WHEN met=1 THEN 0 ELSE 2*canc END
  + 3*msld + 6*met + 6*aids AS charlson_ex_dm,
  diab, diab_c
FROM flag
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE vent24 AS
SELECT c.stay_id,
  MAX(CASE WHEN v.ventilation_status='InvasiveVent'
            AND v.starttime < c.intime + INTERVAL 24 HOUR
            AND v.endtime > c.intime
           THEN 1 ELSE 0 END) AS vent24
FROM first_icu c
LEFT JOIN mimic.mimiciv_derived.ventilation v ON v.stay_id=c.stay_id
GROUP BY c.stay_id
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE vaso24 AS
SELECT c.stay_id,
  MAX(CASE WHEN ie.itemid IN (221906,221289,221662,221653,221749,222315)
            AND ie.starttime >= c.intime
            AND ie.starttime < c.intime + INTERVAL 24 HOUR
           THEN 1 ELSE 0 END) AS vaso24
FROM first_icu c
LEFT JOIN wh.mimic_iv.inputevents ie ON ie.stay_id=c.stay_id
GROUP BY c.stay_id
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE mg24 AS
SELECT c.stay_id,
  MIN(le.valuenum) FILTER (WHERE le.valuenum IS NOT NULL AND le.valuenum>0 AND le.valuenum<10) AS mg_min,
  MAX(le.valuenum) FILTER (WHERE le.valuenum IS NOT NULL AND le.valuenum>0 AND le.valuenum<10) AS mg_max
FROM first_icu c
LEFT JOIN wh.mimic_iv.labevents le
  ON le.subject_id=c.subject_id
 AND le.itemid=50960
 AND le.charttime >= c.intime
 AND le.charttime < c.intime + INTERVAL 24 HOUR
GROUP BY c.stay_id
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE a1c AS
SELECT stay_id, hba1c, hba1c_ge65
FROM (
  SELECT c.stay_id, le.valuenum AS hba1c,
         CASE WHEN le.valuenum>=6.5 THEN 1 ELSE 0 END AS hba1c_ge65,
         ROW_NUMBER() OVER (PARTITION BY c.stay_id
           ORDER BY abs(date_diff('hour', c.admittime, le.charttime)), le.charttime, le.labevent_id, le.valuenum) AS rn
  FROM first_icu c
  INNER JOIN wh.mimic_iv.labevents le ON le.subject_id=c.subject_id AND le.itemid=50852
  WHERE le.valuenum IS NOT NULL AND le.valuenum BETWEEN 3 AND 20
    AND le.charttime >= c.admittime - INTERVAL 90 DAY
    AND le.charttime <= c.admittime + INTERVAL 24 HOUR
)
WHERE rn=1
")

# glucose measurements in first 24h (minute-dedup later in R or SQL)
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE gluc24_raw AS
SELECT c.stay_id, date_trunc('minute', le.charttime) AS tmin, le.valuenum
FROM first_icu c
INNER JOIN wh.mimic_iv.labevents le ON le.subject_id=c.subject_id
WHERE le.itemid IN (50931,50809,52569,52027)
  AND le.valuenum IS NOT NULL AND le.valuenum BETWEEN 10 AND 1000
  AND le.charttime >= c.intime AND le.charttime < c.intime + INTERVAL 24 HOUR
UNION ALL
SELECT c.stay_id, date_trunc('minute', ce.charttime) AS tmin, ce.valuenum
FROM first_icu c
INNER JOIN wh.mimic_iv.chartevents ce ON ce.stay_id=c.stay_id
WHERE ce.itemid IN (220621,226537,225664,228388)
  AND ce.valuenum IS NOT NULL AND ce.valuenum BETWEEN 10 AND 1000
  AND ce.charttime >= c.intime AND ce.charttime < c.intime + INTERVAL 24 HOUR
")
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE gluc24 AS
WITH d AS (
  SELECT stay_id, tmin, avg(valuenum) AS g
  FROM gluc24_raw
  GROUP BY stay_id, tmin
)
SELECT stay_id,
  COUNT(*) AS n_gluc24,
  avg(g) AS gluc_mean,
  stddev_pop(g) AS gluc_sd,
  MIN(g) AS gluc_min,
  MAX(g) AS gluc_max,
  MAX(CASE WHEN g<70 THEN 1 ELSE 0 END) AS hypo70,
  CASE WHEN avg(g)>0 AND stddev_pop(g) IS NOT NULL THEN stddev_pop(g)/avg(g) ELSE NULL END AS gluc_cv
FROM d
GROUP BY stay_id
")

# last rhythm in final 6h before ICU outtime (among survivors used later)
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE last_rhythm AS
SELECT stay_id, last_value, last_is_af
FROM (
  SELECT c.stay_id, e.value AS last_value, e.is_af AS last_is_af,
         ROW_NUMBER() OVER (PARTITION BY c.stay_id ORDER BY e.charttime DESC, e.storetime DESC, e.value, e.valuenum, e.caregiver_id) AS rn
  FROM first_icu c
  INNER JOIN rhythm e ON e.stay_id=c.stay_id
  WHERE e.charttime >= c.outtime - INTERVAL 6 HOUR
    AND e.charttime <= c.outtime
)
WHERE rn=1
")

# early diabetes meds (proxy, not pre-admission) — name match in prescriptions
dbExecute(con, "
CREATE OR REPLACE TEMP TABLE dm_meds AS
SELECT c.stay_id,
  MAX(CASE WHEN regexp_matches(lower(coalesce(pr.drug,'')),
        'insulin|metformin|glipizide|glyburide|glimepiride|sitagliptin|linagliptin|empagliflozin|dapagliflozin|canagliflozin|ertugliflozin|liraglutide|semaglutide|dulaglutide|exenatide|lixisenatide|pioglitazone|rosiglitazone')
      THEN 1 ELSE 0 END) AS dm_meds_early
FROM first_icu c
LEFT JOIN wh.mimic_iv.prescriptions pr
  ON pr.hadm_id=c.hadm_id
 AND pr.starttime <= c.intime + INTERVAL 24 HOUR
 AND (pr.stoptime IS NULL OR pr.stoptime >= c.admittime)
GROUP BY c.stay_id
")

dbExecute(con, "
CREATE OR REPLACE TEMP TABLE analytic AS
SELECT
  x.*,
  CASE WHEN x.has_dm_icd=1 THEN 1 ELSE 0 END AS has_dm,
  CASE WHEN x.has_dm_icd=1 OR coalesce(a.hba1c_ge65,0)=1 OR coalesce(m.dm_meds_early,0)=1 THEN 1 ELSE 0 END AS has_dm_ext,
  a.hba1c, coalesce(a.hba1c_ge65,0) AS hba1c_ge65,
  coalesce(m.dm_meds_early,0) AS dm_meds_early,
  coalesce(ca.csurg_service,0) AS csurg_service,
  CASE WHEN coalesce(ca.csurg_service,0)=1 OR coalesce(ca.cardiac_proc9,0)=1 OR coalesce(ca.cardiac_proc10,0)=1 THEN 1 ELSE 0 END AS cardiac_surgery,
  coalesce(ch.charlson,0) AS charlson,
  coalesce(ch.charlson_ex_dm,0) AS charlson_ex_dm,
  s.sofa, sa.sapsii,
  CASE WHEN sp.stay_id IS NOT NULL THEN 1 ELSE 0 END AS sepsis3,
  coalesce(v.vent24,0) AS vent24,
  coalesce(va.vaso24,0) AS vaso24,
  lab.potassium_min, lab.potassium_max, lab.glucose_min AS lab_glucose_min, lab.glucose_max AS lab_glucose_max,
  mg.mg_min, mg.mg_max,
  g.n_gluc24, g.gluc_mean, g.gluc_sd, g.gluc_min, g.gluc_max, g.gluc_cv, coalesce(g.hypo70,0) AS hypo70,
  lr.last_is_af,
  CASE
    WHEN lower(x.race) LIKE '%white%' AND lower(x.race) NOT LIKE '%black%' THEN 'White'
    WHEN lower(x.race) LIKE '%black%' OR lower(x.race) LIKE '%african%' THEN 'Black'
    WHEN lower(x.race) LIKE '%hispanic%' OR lower(x.race) LIKE '%south american%' THEN 'Hispanic'
    WHEN lower(x.race) LIKE '%asian%' THEN 'Asian'
    ELSE 'Other or Unknown'
  END AS race_coll,
  CASE
    WHEN x.admission_type IN ('EW EMER.','DIRECT EMER.','URGENT') THEN 'Emergency'
    WHEN x.admission_type IN ('ELECTIVE','SURGICAL SAME DAY ADMISSION') THEN 'Elective'
    ELSE 'Other'
  END AS adm_coll
FROM classed2 x
LEFT JOIN a1c a USING (stay_id)
LEFT JOIN dm_meds m USING (stay_id)
LEFT JOIN cardiac ca USING (stay_id)
LEFT JOIN charlson ch USING (stay_id)
LEFT JOIN mimic.mimiciv_derived.first_day_sofa s ON s.stay_id=x.stay_id
LEFT JOIN mimic.mimiciv_derived.sapsii sa ON sa.stay_id=x.stay_id
LEFT JOIN mimic.mimiciv_derived.sepsis3 sp ON sp.stay_id=x.stay_id
LEFT JOIN vent24 v USING (stay_id)
LEFT JOIN vaso24 va USING (stay_id)
LEFT JOIN mimic.mimiciv_derived.first_day_lab lab ON lab.stay_id=x.stay_id
LEFT JOIN mg24 mg USING (stay_id)
LEFT JOIN gluc24 g USING (stay_id)
LEFT JOIN last_rhythm lr USING (stay_id)
")
log("analytic table", dbGetQuery(con,"SELECT COUNT(*) n FROM analytic")$n)

# SHR in subset with both hba1c and a glucose (use first-day mean)
# SHR = glucose_mgdl / (28.7*HbA1c - 46.7)

# ---------------- Pull stay-level frame into memory only (never write, never print) ----------------
log("collect stay-level analytic frame in memory")
dat <- dbGetQuery(con, "
SELECT
  stay_id, subject_id, hadm_id, intime, outtime, dischtime, deathtime,
  anchor_age, gender, anchor_year_group, year_primary, race_coll, adm_coll,
  has_dm_icd, has_dm, has_dm_ext, hba1c, hba1c_ge65, dm_meds_early,
  af_group, af_group_noecg,
  is_preexisting, is_noaf, is_no_af, is_grayzone, is_noaf_broad, is_noaf_2consec, is_noaf_1h,
  t_af_days, n_af, index_af_icd, prior_af_icd, af_pre, early_af_2h,
  af_icu_ecg, has_ecg_icu, af_hosp_ecg, has_ecg_hosp,
  hosp_death, death_28d, death_1y, has_dod, days_to_death, days_to_discharge, icu_los_days,
  CASE WHEN hospital_expire_flag=1 AND outtime >= dischtime - INTERVAL 2 HOUR THEN 1 ELSE 0 END AS died_in_icu,
  has_hf, has_htn, has_ckd, has_stroke, cardiac_surgery, charlson, charlson_ex_dm,
  sofa, sapsii, sepsis3, vent24, vaso24,
  potassium_min, mg_min, n_gluc24, gluc_mean, gluc_cv, hypo70, last_is_af
FROM analytic
")
# derived times (days), floor at 1/24 to keep Surv valid
safe_time <- function(x, default=28) {
  z <- ifelse(is.na(x), default, x)
  pmax(z, 1/24)
}
dat$time_hosp <- safe_time(dat$days_to_discharge, 1)
# if died in hospital, time is days to death/discharge
dat$status_hosp <- as.integer(dat$hosp_death==1)
dat$time_28 <- safe_time(ifelse(dat$death_28d==1, dat$days_to_death, 28), 28)
dat$time_28 <- pmin(dat$time_28, 28)
dat$status_28 <- as.integer(dat$death_28d==1)
dat$time_1y <- safe_time(ifelse(dat$death_1y==1, dat$days_to_death, 365), 365)
dat$time_1y <- pmin(dat$time_1y, 365)
dat$status_1y <- as.integer(dat$death_1y==1)
dat$t_af_days <- ifelse(!is.na(dat$t_af_days), pmax(dat$t_af_days, 1/48), NA_real_)
dat$shr <- ifelse(!is.na(dat$hba1c) & !is.na(dat$gluc_mean) & (28.7*dat$hba1c-46.7)>0,
                  dat$gluc_mean / (28.7*dat$hba1c - 46.7), NA_real_)
dat$af_group <- factor(dat$af_group, levels=c("no_AF","NOAF","preexisting","grayzone","indeterminate"))
dat$sex <- factor(dat$gender, levels=c("F","M"))
dat$race_coll <- factor(dat$race_coll, levels=c("White","Black","Hispanic","Asian","Other or Unknown"))
dat$adm_coll <- factor(dat$adm_coll, levels=c("Emergency","Elective","Other"))
dat$year_f <- factor(dat$anchor_year_group)

# primary analysis subset flags (do not print dat)
in3 <- dat$af_group %in% c("no_AF","NOAF","preexisting")
in3p <- in3 & dat$year_primary==1

# ---------------- Flow ----------------
log("STROBE flow counts")
flow <- dbGetQuery(con, "
SELECT * FROM (
  SELECT '01_adult_first_icu' AS step, COUNT(*) AS n, COUNT(*) FILTER (WHERE has_dm_icd=1) AS n_dm FROM first_icu
  UNION ALL
  SELECT '02_year_2008_2019', COUNT(*) FILTER (WHERE year_primary=1), COUNT(*) FILTER (WHERE year_primary=1 AND has_dm_icd=1) FROM first_icu
  UNION ALL
  SELECT '03_year_2020_2022', COUNT(*) FILTER (WHERE year_primary=0), COUNT(*) FILTER (WHERE year_primary=0 AND has_dm_icd=1) FROM first_icu
  UNION ALL
  SELECT '04_has_rhythm_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND has_rhythm=1), COUNT(*) FILTER (WHERE year_primary=1 AND has_rhythm=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '05_preexisting_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '06_NOAF_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_noaf=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_noaf=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '07_no_AF_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_no_af=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_no_af=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '08_grayzone_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_grayzone=1), COUNT(*) FILTER (WHERE year_primary=1 AND is_grayzone=1 AND has_dm_icd=1) FROM classed2
  UNION ALL
  SELECT '09_indeterminate_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting=0 AND is_noaf=0 AND is_no_af=0 AND is_grayzone=0),
         COUNT(*) FILTER (WHERE year_primary=1 AND has_dm_icd=1 AND is_preexisting=0 AND is_noaf=0 AND is_no_af=0 AND is_grayzone=0) FROM classed2
  UNION ALL
  SELECT '10_three_group_primary_2008_2019', COUNT(*) FILTER (WHERE year_primary=1 AND is_preexisting+is_noaf+is_no_af=1),
         COUNT(*) FILTER (WHERE year_primary=1 AND has_dm_icd=1 AND is_preexisting+is_noaf+is_no_af=1) FROM classed2
)
ORDER BY step
")
write_both(flow, "01_flow_counts.csv", c("n","n_dm"))

exp_n <- dbGetQuery(con, "
SELECT year_primary, has_dm_icd AS has_dm, af_group, COUNT(*) AS n,
       COUNT(*) FILTER (WHERE hosp_death=1) AS n_hosp_death,
       COUNT(*) FILTER (WHERE death_28d=1) AS n_death_28d,
       COUNT(*) FILTER (WHERE death_1y=1) AS n_death_1y
FROM classed2
GROUP BY 1,2,3
ORDER BY 1 DESC, 2 DESC, 3
")
write_both(exp_n, "02_group_counts_mortality.csv", c("n","n_hosp_death","n_death_28d","n_death_1y"))

# ---------------- Missingness ----------------
log("missingness table")
miss_vars <- c("sofa","sapsii","potassium_min","mg_min","gluc_mean","hba1c","race_coll")
miss_rows <- list()
sub <- dat[in3p, ]
for (g in c("no_AF","NOAF","preexisting")) {
  for (dm in c(0,1)) {
    s <- sub[sub$af_group==g & sub$has_dm==dm, ]
    row <- data.frame(af_group=g, has_dm=dm, n=nrow(s), stringsAsFactors=FALSE)
    for (v in c("sofa","sapsii","potassium_min","mg_min","gluc_mean","hba1c")) {
      row[[paste0("n_miss_", v)]] <- sum(is.na(s[[v]]))
      row[[paste0("pct_miss_", v)]] <- if (nrow(s)>0) round(100*mean(is.na(s[[v]])), 1) else NA
    }
    miss_rows[[length(miss_rows)+1]] <- row
  }
}
miss <- do.call(rbind, miss_rows)
write_both(miss, "03_missingness_by_group.csv", grep("^n", names(miss), value=TRUE))

# ---------------- Table 1 ----------------
log("Table 1")
tab1_src <- dat[in3p, ]
tab1_src$af_dm <- paste0(ifelse(tab1_src$has_dm==1,"DM: ","nonDM: "), tab1_src$af_group)
# variables
tab1_src$gluc_mean_r <- tab1_src$gluc_mean
if (requireNamespace("tableone", quietly=TRUE)) {
  vars <- c("anchor_age","sex","race_coll","adm_coll","anchor_year_group",
            "sofa","sapsii","charlson_ex_dm","sepsis3","has_hf","has_htn","has_ckd",
            "cardiac_surgery","vent24","vaso24","potassium_min","mg_min","gluc_mean","hypo70",
            "hosp_death","death_28d","death_1y","icu_los_days")
  factorVars <- c("sex","race_coll","adm_coll","anchor_year_group","sepsis3","has_hf","has_htn","has_ckd",
                  "cardiac_surgery","vent24","vaso24","hypo70","hosp_death","death_28d","death_1y")
  t1 <- tableone::CreateTableOne(vars=vars, strata="af_dm", data=tab1_src, factorVars=factorVars, addOverall=FALSE, test=TRUE)
  t1m <- print(t1, printToggle=FALSE, showAllLevels=TRUE, quote=FALSE, noSpaces=TRUE, smd=TRUE)
  t1df <- data.frame(item=rownames(t1m), t1m, check.names=FALSE, stringsAsFactors=FALSE)
  # TableOne print can expose cells <=9; keep raw print INTERNAL only. Public uses 04_table1_counts.csv (masked).
  write_internal(t1df, "04_table1_raw.csv")
} else {
  log("tableone not available unexpectedly")
}

# compact numeric Table 1 for masking control
t1_counts <- do.call(rbind, lapply(split(tab1_src, list(tab1_src$has_dm, tab1_src$af_group), drop=TRUE), function(s) {
  data.frame(
    has_dm=s$has_dm[1], af_group=as.character(s$af_group[1]), n=nrow(s),
    age_median=median(s$anchor_age, na.rm=TRUE),
    age_iqr_l=as.numeric(quantile(s$anchor_age, 0.25, na.rm=TRUE)),
    age_iqr_h=as.numeric(quantile(s$anchor_age, 0.75, na.rm=TRUE)),
    n_female=sum(s$gender=="F", na.rm=TRUE),
    sofa_median=median(s$sofa, na.rm=TRUE),
    saps_median=median(s$sapsii, na.rm=TRUE),
    n_sepsis=sum(s$sepsis3==1, na.rm=TRUE),
    n_hf=sum(s$has_hf==1), n_htn=sum(s$has_htn==1), n_ckd=sum(s$has_ckd==1),
    n_cardiac=sum(s$cardiac_surgery==1), n_vent=sum(s$vent24==1), n_vaso=sum(s$vaso24==1),
    n_hypo=sum(s$hypo70==1, na.rm=TRUE),
    n_hosp_death=sum(s$hosp_death==1), n_death_28d=sum(s$death_28d==1), n_death_1y=sum(s$death_1y==1),
    icu_los_median=median(s$icu_los_days, na.rm=TRUE),
    stringsAsFactors=FALSE)
}))
t1_counts$pct_hosp_death <- ifelse(t1_counts$n>0, round(100*t1_counts$n_hosp_death/t1_counts$n, 1), NA)
t1_counts$pct_death_28d <- ifelse(t1_counts$n>0, round(100*t1_counts$n_death_28d/t1_counts$n, 1), NA)
t1_counts$pct_death_1y <- ifelse(t1_counts$n>0, round(100*t1_counts$n_death_1y/t1_counts$n, 1), NA)
write_both(t1_counts, "04_table1_counts.csv", grep("^n", names(t1_counts), value=TRUE))

# NOAF incidence among at-risk (not preexisting), 2008-2019
inc <- dbGetQuery(con, "
SELECT has_dm_icd AS has_dm,
  COUNT(*) FILTER (WHERE is_preexisting=0 AND (is_noaf=1 OR is_no_af=1 OR is_grayzone=1 OR (has_rhythm=1))) AS n_atrisk_broad,
  COUNT(*) FILTER (WHERE is_noaf=1) AS n_noaf,
  COUNT(*) FILTER (WHERE is_preexisting=1) AS n_pre,
  COUNT(*) FILTER (WHERE year_primary=1) AS n_year
FROM classed2
WHERE year_primary=1
GROUP BY 1
")
# cleaner: at-risk = not preexisting and classifiable as NOAF or no_AF (exclude gray/indet from denominator for incidence among 'at risk for new AF')
inc2 <- with(dat[dat$year_primary==1,], {
  do.call(rbind, lapply(c(0,1), function(dm) {
    s <- dat[dat$year_primary==1 & dat$has_dm==dm, ]
    n_pre <- sum(s$af_group=="preexisting")
    n_noaf <- sum(s$af_group=="NOAF")
    n_no <- sum(s$af_group=="no_AF")
    n_atrisk <- n_noaf + n_no
    data.frame(has_dm=dm, n_total=nrow(s), n_preexisting=n_pre, n_NOAF=n_noaf, n_no_AF=n_no,
               n_atrisk=n_atrisk,
               inc_pct_atrisk=ifelse(n_atrisk>0, round(100*n_noaf/n_atrisk, 1), NA),
               inc_pct_all=round(100*n_noaf/nrow(s), 1),
               stringsAsFactors=FALSE)
  }))
})
write_both(inc2, "05_noaf_incidence.csv", c("n_total","n_preexisting","n_NOAF","n_no_AF","n_atrisk"))
log("flow, missing, table1, incidence written")

# ---------------- Descriptive KM (not the primary causal estimator) ----------------
log("Kaplan-Meier figures")
km_dat <- dat[in3p, ]
km_dat$grp <- factor(km_dat$af_group, levels=c("no_AF","NOAF","preexisting"),
                     labels=c("No AF","NOAF","Pre-existing AF"))
km_dat$dm_lab <- factor(km_dat$has_dm, levels=c(0,1), labels=c("Without diabetes","With diabetes"))
if (requireNamespace("survminer", quietly=TRUE)) {
  make_km <- function(time, status, ylab, fname) {
    fit <- survfit(Surv(km_dat[[time]], km_dat[[status]]) ~ grp + dm_lab, data=km_dat)
    p <- survminer::ggsurvplot(
      fit, data=km_dat, conf.int=TRUE, risk.table=TRUE, risk.table.height=0.28,
      palette=c(okabe["no_AF"], okabe["NOAF"], okabe["preexisting"],
                okabe["no_AF"], okabe["NOAF"], okabe["preexisting"]),
      linetype="solid",
      legend.title="", xlab="Days from ICU admission", ylab=ylab,
      ggtheme=theme_journal(11),
      fontsize=3.5, tables.theme=survminer::theme_cleantable()
    )
    # survminer returns a list; arrange and save
    p$plot <- p$plot + facet_wrap(~dm_lab, ncol=2) + theme_journal(11) +
      labs(caption="Descriptive KM by baseline group. Primary inference uses time-varying Cox (SAP v1).")
    # Simpler custom KM by diabetes facet to control colors per panel
    p2 <- ggplot(NULL) # placeholder replaced below
    invisible(fit)
  }
}
# Custom faceted KM for control of style
plot_km <- function(df, time, status, title, fname) {
  d <- df
  d$time <- d[[time]]; d$status <- d[[status]]
  d <- d[!is.na(d$time) & !is.na(d$status) & !is.na(d$grp), ]
  fit <- survfit(Surv(time, status) ~ grp + dm_lab, data=d)
  s <- survminer::surv_summary(fit, data=d)
  # strata labels
  p <- ggplot(s, aes(x=time, y=surv, color=grp, fill=grp)) +
    geom_step(linewidth=0.6) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.15, color=NA) +
    facet_wrap(~dm_lab, ncol=2) +
    scale_color_manual(values=c("No AF"=okabe[["no_AF"]], "NOAF"=okabe[["NOAF"]], "Pre-existing AF"=okabe[["preexisting"]])) +
    scale_fill_manual(values=c("No AF"=okabe[["no_AF"]], "NOAF"=okabe[["NOAF"]], "Pre-existing AF"=okabe[["preexisting"]])) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
    labs(x="Days from ICU admission", y="Survival probability", color=NULL, fill=NULL,
         title=title,
         caption="Descriptive only; groups are baseline labels. Primary model is time-varying Cox (SAP v1.0).") +
    theme_journal(11)
  save_fig(p, fname, width=8.5, height=5.2)
  p
}
plot_km(km_dat, "time_28", "status_28", "28-day all-cause survival (descriptive)", "fig2_km_28d")
plot_km(km_dat, "time_1y", "status_1y", "1-year all-cause survival (descriptive)", "fig2_km_1y")
plot_km(km_dat, "time_hosp", "status_hosp", "In-hospital survival (descriptive; censored at discharge)", "fig2_km_hosp")

# ---------------- Time-varying long format ----------------
# Build counting-process data without writing to disk
make_tvc <- function(df, time_var, status_var) {
  d <- df
  d$tstop0 <- d[[time_var]]
  d$ev <- d[[status_var]]
  d <- d[!is.na(d$tstop0) & d$tstop0>0, ]
  d$g <- as.character(d$af_group)
  split_noaf <- d$g=="NOAF" & !is.na(d$t_af_days) & d$t_af_days < d$tstop0
  a <- d
  a$tstart <- 0
  a$tstop <- ifelse(split_noaf, d$t_af_days, d$tstop0)
  a$event <- ifelse(split_noaf, 0L, d$ev)
  a$exp <- ifelse(d$g=="NOAF", "no_AF", d$g)
  b <- d[split_noaf, , drop=FALSE]
  if (nrow(b)) {
    b$tstart <- b$t_af_days
    b$tstop <- b$tstop0
    b$event <- b$ev
    b$exp <- "NOAF"
    out <- rbind(a, b)
  } else out <- a
  out$exp <- factor(out$exp, levels=c("no_AF","NOAF","preexisting"))
  out[out$tstop > out$tstart, ]
}

covars <- c("anchor_age","sex","race_coll","adm_coll","year_f","sofa","sapsii",
            "charlson_ex_dm","sepsis3","has_hf","has_htn","has_ckd","cardiac_surgery",
            "vent24","vaso24","potassium_min","mg_min","gluc_mean")
form_rhs <- paste(c("exp", covars), collapse=" + ")

extract_hr <- function(fit, model, outcome, subset_lab, n=NULL, n_event=NULL) {
  sm <- summary(fit)$conf.int
  if (is.null(sm)) return(NULL)
  rn <- rownames(sm)
  keep <- grepl("^exp", rn)
  if (!any(keep)) return(NULL)
  data.frame(
    model=model, outcome=outcome, subset=subset_lab, term=sub("^exp","", rn[keep]),
    hr=sm[keep, "exp(coef)"],
    lcl=sm[keep, "lower .95"],
    ucl=sm[keep, "upper .95"],
    p=summary(fit)$coefficients[keep, "Pr(>|z|)"],
    n=if (is.null(n)) fit$n else n,
    n_event=if (is.null(n_event)) as.integer(fit$nevent) else as.integer(n_event),
    stringsAsFactors=FALSE
  )
}

# t_af_days is missing by design in no_AF; requiring it in complete.cases drops the reference group.
drop_noaf_missing_onset <- function(d, lab="") {
  noaf_miss <- !is.na(d$af_group) & as.character(d$af_group)=="NOAF" & (is.na(d$t_af_days) | d$t_af_days<=0)
  if (any(noaf_miss)) {
    log("drop NOAF missing t_af_days:", sum(noaf_miss), lab)
    d <- d[!noaf_miss, ]
  }
  d
}

fit_tvc_cox <- function(df, time_var, status_var, outcome, subset_lab, use_imp=FALSE) {
  d <- df
  # complete covariates for complete-case; do not require t_af_days for all groups
  need <- c(time_var, status_var, "af_group", covars)
  if (!use_imp) {
    d <- d[stats::complete.cases(d[, intersect(need, names(d)), drop=FALSE]), ]
  }
  d <- drop_noaf_missing_onset(d, paste(outcome, subset_lab))
  if (nrow(d)<50 || sum(d[[status_var]], na.rm=TRUE)<10) {
    log("skip cox (small):", outcome, subset_lab, "n", nrow(d))
    return(NULL)
  }
  long <- make_tvc(d, time_var, status_var)
  f <- as.formula(paste("Surv(tstart, tstop, event) ~", form_rhs))
  fit <- try(coxph(f, data=long, ties="efron"), silent=TRUE)
  if (inherits(fit, "try-error")) { log("cox failed", outcome, subset_lab, as.character(fit)); return(NULL) }
  extract_hr(fit, "tvc_complete_case", outcome, subset_lab,
             n=nrow(d), n_event=sum(d[[status_var]], na.rm=TRUE))
}

# Primary complete-case TVC (also run mice below)
log("time-varying Cox complete-case (primary years, three groups)")
hr_list <- list()
for (dm_lab in list(list(lab="DM", fun=function(d) d$has_dm==1),
                    list(lab="nonDM", fun=function(d) d$has_dm==0),
                    list(lab="all", fun=function(d) TRUE))) {
  dsub <- dat[in3p & dm_lab$fun(dat), ]
  for (oc in list(c("time_hosp","status_hosp","hosp_death"),
                  c("time_28","status_28","death_28d"),
                  c("time_1y","status_1y","death_1y"))) {
    hr_list[[length(hr_list)+1]] <- fit_tvc_cox(dsub, oc[1], oc[2], oc[3], dm_lab$lab, use_imp=FALSE)
  }
}
# Logistic hospital death complete-case
log("logistic hospital death")
for (dm_lab in c("DM","nonDM","all")) {
  dsub <- dat[in3p, ]
  if (dm_lab=="DM") dsub <- dsub[dsub$has_dm==1, ]
  if (dm_lab=="nonDM") dsub <- dsub[dsub$has_dm==0, ]
  dsub$exp <- factor(dsub$af_group, levels=c("no_AF","NOAF","preexisting"))
  dcc <- dsub[stats::complete.cases(dsub[, c("hosp_death","exp", covars)]), ]
  if (nrow(dcc)<50) next
  fit <- try(glm(as.formula(paste("hosp_death ~ exp +", paste(covars, collapse=" + "))),
                 data=dcc, family=binomial()), silent=TRUE)
  if (inherits(fit, "try-error")) next
  sm <- summary(fit)$coefficients
  ci <- try(confint.default(fit), silent=TRUE)
  terms <- grep("^exp", rownames(sm), value=TRUE)
  for (tm in terms) {
    hr_list[[length(hr_list)+1]] <- data.frame(
      model="logistic_complete_case", outcome="hosp_death_OR", subset=dm_lab, term=sub("^exp","",tm),
      hr=exp(sm[tm,"Estimate"]),
      lcl=exp(ci[tm,1]), ucl=exp(ci[tm,2]),
      p=sm[tm,"Pr(>|z|)"], n=nrow(dcc), n_event=sum(dcc$hosp_death),
      stringsAsFactors=FALSE)
  }
}

# ---------------- mice imputation then TVC (primary adjusted model) ----------------
log("mice m=10 on baseline covariates, then tvc cox (Rubin pool)")
imp_vars <- c("sofa","sapsii","potassium_min","mg_min","gluc_mean")
# race/adm/year/sex complete as factors (Unknown used)
do_mice <- TRUE
if (!requireNamespace("mice", quietly=TRUE)) {
  log("mice not installed; skip MI (should not happen)")
  do_mice <- FALSE
}
hr_mi <- list()
if (do_mice) {
  dmi0 <- dat[in3p, c("stay_id","has_dm","af_group","t_af_days",
                      "time_hosp","status_hosp","time_28","status_28","time_1y","status_1y",
                      covars)]
  # mice on covariates only + predictors
  pred_df <- dmi0[, c("has_dm","af_group", covars)]
  pred_df$sepsis3 <- factor(pred_df$sepsis3)
  pred_df$has_hf <- factor(pred_df$has_hf)
  pred_df$has_htn <- factor(pred_df$has_htn)
  pred_df$has_ckd <- factor(pred_df$has_ckd)
  pred_df$cardiac_surgery <- factor(pred_df$cardiac_surgery)
  pred_df$vent24 <- factor(pred_df$vent24)
  pred_df$vaso24 <- factor(pred_df$vaso24)
  ini <- try(mice::mice(pred_df, maxit=0, print=FALSE), silent=TRUE)
  if (inherits(ini, "try-error")) {
    log("mice setup failed", as.character(ini)); do_mice <- FALSE
  } else {
    meth <- ini$method
    # only impute continuous labs/severity; leave complete factors
    for (v in names(meth)) meth[v] <- ""
    for (v in imp_vars) if (v %in% names(meth)) meth[v] <- "pmm"
    log("running mice ...")
    imp <- try(mice::mice(pred_df, m=10, maxit=8, method=meth, seed=SEED, print=FALSE), silent=TRUE)
    if (inherits(imp, "try-error")) {
      log("mice failed", as.character(imp)); do_mice <- FALSE
    } else {
      for (oc in list(c("time_28","status_28","death_28d"),
                      c("time_hosp","status_hosp","hosp_death"),
                      c("time_1y","status_1y","death_1y"))) {
        for (dm_lab in c("DM","nonDM","all")) {
          coefs <- list(); vcovs <- list(); ns <- c(); evs <- c()
          for (k in 1:10) {
            comp <- mice::complete(imp, k)
            dfill <- dmi0
            for (v in imp_vars) dfill[[v]] <- comp[[v]]
            if (dm_lab=="DM") dfill <- dfill[dfill$has_dm==1, ]
            if (dm_lab=="nonDM") dfill <- dfill[dfill$has_dm==0, ]
            long <- make_tvc(dfill, oc[1], oc[2])
            # restore factor types
            long$sex <- factor(long$sex)
            long$race_coll <- factor(long$race_coll)
            long$adm_coll <- factor(long$adm_coll)
            long$year_f <- factor(long$year_f)
            f <- as.formula(paste("Surv(tstart, tstop, event) ~", form_rhs))
            fit <- try(coxph(f, data=long, ties="efron"), silent=TRUE)
            if (inherits(fit, "try-error")) next
            cf <- coef(fit); vc <- vcov(fit)
            keep <- grepl("^exp", names(cf))
            coefs[[length(coefs)+1]] <- cf[keep]
            vcovs[[length(vcovs)+1]] <- vc[keep, keep, drop=FALSE]
            ns <- c(ns, nrow(dfill)); evs <- c(evs, sum(dfill[[oc[2]]], na.rm=TRUE))
          }
          if (length(coefs)<3) { log("too few MI fits", oc[3], dm_lab); next }
          # Rubin pool for each exp term
          tn <- names(coefs[[1]])
          for (tm in tn) {
            qs <- sapply(coefs, function(b) unname(b[tm]))
            ubar <- mean(sapply(seq_along(vcovs), function(i) vcovs[[i]][tm,tm]))
            bvar <- var(qs)
            m <- length(qs)
            Tvar <- ubar + (1+1/m)*bvar
            se <- sqrt(Tvar)
            est <- mean(qs)
            z <- est/se
            p <- 2*pnorm(-abs(z))
            hr_mi[[length(hr_mi)+1]] <- data.frame(
              model="tvc_mice_m10", outcome=oc[3], subset=dm_lab, term=sub("^exp","",tm),
              hr=exp(est), lcl=exp(est-1.96*se), ucl=exp(est+1.96*se), p=p,
              n=round(mean(ns)), n_event=round(mean(evs)),
              stringsAsFactors=FALSE)
          }
        }
      }
      log("mice tvc pooled")
    }
  }
}

hr_all <- do.call(rbind, c(hr_list[!sapply(hr_list, is.null)], hr_mi))
if (!is.null(hr_all) && nrow(hr_all)) {
  hr_all$evalue <- mapply(evalue_hr, hr_all$hr)
  hr_all$evalue_lcl <- mapply(function(hr,lcl,ucl) {
    # E-value for the CI bound closer to 1
    if (is.na(hr)||is.na(lcl)||is.na(ucl)) return(NA)
    if (hr>=1) evalue_hr(lcl) else evalue_hr(ucl)
  }, hr_all$hr, hr_all$lcl, hr_all$ucl)
  write_both(hr_all, "06_hr_main_and_logistic.csv", c("n","n_event"))
  write_internal(hr_all, "06_hr_main_and_logistic.csv")
  log("HR table rows", nrow(hr_all))
} else {
  log("WARNING no HR rows")
}

# ---------------- Landmark 48h ----------------
log("landmark 48h")
lm_dat <- dat[in3p & dat$time_28 > 2, ]  # alive (or still followed) beyond 2 days; exclude early deaths
# For 28d follow-up starting at day 2: those who died before 48h are excluded
lm_dat <- lm_dat[!(lm_dat$status_28==1 & lm_dat$time_28<=2), ]
lm_dat$lm_group <- as.character(lm_dat$af_group)
# NOAF only if t_af <= 2 days; if t_af>2, treat as no_AF at landmark
lm_dat$lm_group[lm_dat$af_group=="NOAF" & (is.na(lm_dat$t_af_days) | lm_dat$t_af_days>2)] <- "no_AF"
lm_dat$lm_group <- factor(lm_dat$lm_group, levels=c("no_AF","NOAF","preexisting"))
lm_dat$lm_time28 <- pmax(lm_dat$time_28 - 2, 1/24)
lm_dat$lm_time1y <- pmax(lm_dat$time_1y - 2, 1/24)
lm_n <- aggregate(stay_id~has_dm+lm_group, lm_dat, length)
names(lm_n)[3] <- "n"
write_both(lm_n, "07_landmark48_counts.csv", "n")

fit_lm <- function(df, time, status, outcome, subset) {
  d <- df[stats::complete.cases(df[, c(time,status,"lm_group", covars)]), ]
  d$exp <- d$lm_group
  if (nrow(d)<50) return(NULL)
  f <- as.formula(paste("Surv(", time, ",", status, ") ~ exp +", paste(covars, collapse=" + ")))
  fit <- try(coxph(f, data=d, ties="efron"), silent=TRUE)
  if (inherits(fit,"try-error")) return(NULL)
  extract_hr(fit, "landmark_48h_cc", outcome, subset)
}
hr_lm <- list()
for (dm_lab in c("DM","nonDM","all")) {
  dsub <- lm_dat
  if (dm_lab=="DM") dsub <- dsub[dsub$has_dm==1,]
  if (dm_lab=="nonDM") dsub <- dsub[dsub$has_dm==0,]
  hr_lm[[length(hr_lm)+1]] <- fit_lm(dsub, "lm_time28", "status_28", "death_28d", dm_lab)
  hr_lm[[length(hr_lm)+1]] <- fit_lm(dsub, "lm_time1y", "status_1y", "death_1y", dm_lab)
}
hr_lm <- do.call(rbind, hr_lm[!sapply(hr_lm,is.null)])
if (!is.null(hr_lm) && nrow(hr_lm)) {
  hr_lm$evalue <- mapply(evalue_hr, hr_lm$hr)
  write_both(hr_lm, "07_hr_landmark48.csv", c("n","n_event"))
}

# ---------------- Interaction and RERI (complete-case TVC, exclude preexisting for 2x2 RERI) ----------------
log("interaction and RERI")
int_dat <- dat[in3p & dat$af_group %in% c("no_AF","NOAF"), ]
int_dat <- int_dat[stats::complete.cases(int_dat[, c("time_28","status_28","af_group","has_dm", covars)]), ]
int_dat <- drop_noaf_missing_onset(int_dat, "RERI 28d")
long_i <- make_tvc(int_dat, "time_28", "status_28")
long_i$dm <- long_i$has_dm
f_int <- as.formula(paste("Surv(tstart, tstop, event) ~ exp * dm +", paste(covars, collapse=" + ")))
fit_int <- try(coxph(f_int, data=long_i, ties="efron"), silent=TRUE)
reri_rows <- NULL
if (!inherits(fit_int, "try-error")) {
  sm <- summary(fit_int)
  # terms: expNOAF, dm, expNOAF:dm
  cf <- coef(fit_int); vc <- vcov(fit_int)
  b_e <- unname(cf["expNOAF"]); b_d <- unname(cf["dm"]); b_ed <- unname(cf["expNOAF:dm"])
  if (!any(is.na(c(b_e,b_d,b_ed)))) {
    hr_e <- exp(b_e); hr_d <- exp(b_d); hr_ed <- exp(b_e + b_d + b_ed)
    reri <- hr_ed - hr_e - hr_d + 1
    ap <- reri / hr_ed
    sidx <- hr_ed / (hr_e + hr_d - 1)
    # delta method for RERI (Hosmer)
    # dRERI/db = (hr_ed * x) - (hr_e for e only) etc.
    ge <- hr_ed - hr_e
    gd <- hr_ed - hr_d
    ged <- hr_ed
    g <- c(ge, gd, ged)
    idx <- c("expNOAF","dm","expNOAF:dm")
    V <- vc[idx, idx, drop=FALSE]
    se_reri <- sqrt(as.numeric(t(g) %*% V %*% g))
    reri_rows <- data.frame(
      outcome="death_28d", n=nrow(int_dat), n_event=sum(int_dat$status_28, na.rm=TRUE),
      hr_noaf_nondm=hr_e, hr_dm_no_noaf=hr_d, hr_both=hr_ed,
      reri=reri, reri_lcl=reri-1.96*se_reri, reri_ucl=reri+1.96*se_reri,
      ap=ap, synergy=sidx,
      p_interaction_multiplicative=sm$coefficients["expNOAF:dm","Pr(>|z|)"],
      stringsAsFactors=FALSE
    )
    log("RERI 28d", round(reri,3), "Pint", signif(reri_rows$p_interaction_multiplicative,3))
  }
  # also hospital and 1y
  write_both(extract_hr(fit_int, "tvc_interaction_cc", "death_28d", "NOAF_x_DM",
                         n=nrow(int_dat), n_event=sum(int_dat$status_28, na.rm=TRUE)),
             "08_hr_interaction_terms_28d.csv", c("n","n_event"))
}
if (!is.null(reri_rows)) write_both(reri_rows, "08_reri_28d.csv", c("n","n_event"))

# Repeat RERI for hosp and 1y briefly
reri_all <- list(reri_rows)
for (oc in list(c("time_hosp","status_hosp","hosp_death"), c("time_1y","status_1y","death_1y"))) {
  d2 <- dat[in3p & dat$af_group %in% c("no_AF","NOAF"), ]
  d2 <- d2[stats::complete.cases(d2[, c(oc[1],oc[2],"af_group","has_dm", covars)]), ]
  d2 <- drop_noaf_missing_onset(d2, oc[3])
  long2 <- make_tvc(d2, oc[1], oc[2]); long2$dm <- long2$has_dm
  fit2 <- try(coxph(as.formula(paste("Surv(tstart, tstop, event) ~ exp * dm +", paste(covars, collapse=" + "))),
                    data=long2, ties="efron"), silent=TRUE)
  if (inherits(fit2,"try-error")) next
  cf <- coef(fit2); vc <- vcov(fit2)
  if (!all(c("expNOAF","dm","expNOAF:dm") %in% names(cf))) next
  b_e <- unname(cf["expNOAF"]); b_d <- unname(cf["dm"]); b_ed <- unname(cf["expNOAF:dm"])
  hr_e <- exp(b_e); hr_d <- exp(b_d); hr_ed <- exp(b_e+b_d+b_ed)
  reri <- hr_ed-hr_e-hr_d+1
  ge <- hr_ed-hr_e; gd <- hr_ed-hr_d; ged <- hr_ed
  V <- vc[c("expNOAF","dm","expNOAF:dm"), c("expNOAF","dm","expNOAF:dm")]
  se <- sqrt(as.numeric(t(c(ge,gd,ged)) %*% V %*% c(ge,gd,ged)))
  reri_all[[length(reri_all)+1]] <- data.frame(
    outcome=oc[3], n=nrow(d2), n_event=sum(d2[[oc[2]]], na.rm=TRUE),
    hr_noaf_nondm=hr_e, hr_dm_no_noaf=hr_d, hr_both=hr_ed,
    reri=reri, reri_lcl=reri-1.96*se, reri_ucl=reri+1.96*se,
    ap=reri/hr_ed, synergy=hr_ed/(hr_e+hr_d-1),
    p_interaction_multiplicative=summary(fit2)$coefficients["expNOAF:dm","Pr(>|z|)"],
    stringsAsFactors=FALSE)
}
reri_bind <- do.call(rbind, reri_all[!sapply(reri_all,is.null)])
if (!is.null(reri_bind) && nrow(reri_bind)) write_both(reri_bind, "08_reri_all_outcomes.csv", c("n","n_event"))

# ---------------- Sensitivities (complete-case TVC, 28d, DM and all) ----------------
log("sensitivity models")
sens_hr <- list()

# helper to recode group then fit 28d
fit_alt_group <- function(df, group_col, label, subset_lab="DM") {
  d <- df
  d$af_group <- factor(d[[group_col]], levels=c("no_AF","NOAF","preexisting"))
  d <- d[!is.na(d$af_group), ]
  fit_tvc_cox(d, "time_28", "status_28", "death_28d", paste(subset_lab, label), use_imp=FALSE)
}

# 1) ECG-confirmed NOAF: drop NOAF without in-hospital window ECG AF; keep others
d_ecg <- dat[in3p, ]
d_ecg <- d_ecg[!(d_ecg$af_group=="NOAF" & d_ecg$af_hosp_ecg==0), ]
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ecg[d_ecg$has_dm==1,], "time_28","status_28","death_28d","DM ECG-confirmed NOAF")
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ecg, "time_28","status_28","death_28d","all ECG-confirmed NOAF")

# 2) exclude cardiac surgery
d_nocs <- dat[in3p & dat$cardiac_surgery==0, ]
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_nocs[d_nocs$has_dm==1,], "time_28","status_28","death_28d","DM excl cardiac surgery")
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_nocs, "time_28","status_28","death_28d","all excl cardiac surgery")

# 3) include 2020-2022
d_ally <- dat[in3, ]
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ally[d_ally$has_dm==1,], "time_28","status_28","death_28d","DM incl 2020-2022")
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ally, "time_28","status_28","death_28d","all incl 2020-2022")

# 4) gray zone reclassification
# A: gray -> preexisting
d_ga <- dat[dat$year_primary==1 & dat$af_group %in% c("no_AF","NOAF","preexisting","grayzone"), ]
d_ga$af_group <- factor(ifelse(as.character(d_ga$af_group)=="grayzone","preexisting", as.character(d_ga$af_group)),
                        levels=c("no_AF","NOAF","preexisting"))
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ga[d_ga$has_dm==1,], "time_28","status_28","death_28d","DM gray->preexisting")
# B: gray -> no_AF
d_gb <- dat[dat$year_primary==1 & dat$af_group %in% c("no_AF","NOAF","preexisting","grayzone"), ]
d_gb$af_group <- factor(ifelse(as.character(d_gb$af_group)=="grayzone","no_AF", as.character(d_gb$af_group)),
                        levels=c("no_AF","NOAF","preexisting"))
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_gb[d_gb$has_dm==1,], "time_28","status_28","death_28d","DM gray->no_AF")

# 5) broader NOAF / 2consec / 1h : rebuild three-group among year_primary using alt flags
alt_three <- function(is_noaf_col, lab) {
  d <- dat[dat$year_primary==1, ]
  g <- rep(NA_character_, nrow(d))
  g[d$is_preexisting==1] <- "preexisting"
  g[d[[is_noaf_col]]==1] <- "NOAF"
  # no_AF if not preexisting, not this NOAF, and original no_AF (never AF, not gray)
  g[is.na(g) & d$is_no_af==1] <- "no_AF"
  d$af_group <- factor(g, levels=c("no_AF","NOAF","preexisting"))
  d <- d[!is.na(d$af_group), ]
  fit_tvc_cox(d[d$has_dm==1,], "time_28","status_28","death_28d", paste("DM", lab))
}
sens_hr[[length(sens_hr)+1]] <- alt_three("is_noaf_broad", "broad NOAF")
sens_hr[[length(sens_hr)+1]] <- alt_three("is_noaf_2consec", "NOAF 2consec")
sens_hr[[length(sens_hr)+1]] <- alt_three("is_noaf_1h", "NOAF >=1h")

# 6) definition without pre-ICU ECG (old preexisting)
d_ne <- dat[dat$year_primary==1, ]
d_ne$af_group <- factor(d_ne$af_group_noecg, levels=c("no_AF","NOAF","preexisting"))
d_ne <- d_ne[!is.na(d_ne$af_group), ]
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ne[d_ne$has_dm==1,], "time_28","status_28","death_28d","DM no pre-ICU ECG in def")

# 7) extended diabetes definition: compare using has_dm_ext in 'DM' subset
d_ext <- dat[in3p & dat$has_dm_ext==1, ]
sens_hr[[length(sens_hr)+1]] <- fit_tvc_cox(d_ext, "time_28","status_28","death_28d","DM-extended ICD/A1c/meds")

sens_bind <- do.call(rbind, sens_hr[!sapply(sens_hr,is.null)])
if (!is.null(sens_bind) && nrow(sens_bind)) {
  sens_bind$evalue <- mapply(evalue_hr, sens_bind$hr)
  write_both(sens_bind, "09_hr_sensitivity.csv", c("n","n_event"))
}

# Glucose strata and hypoglycemia within DM, 28d complete-case TVC
log("DM glucose strata")
dm28 <- dat[in3p & dat$has_dm==1, ]
# tertiles of CV and mean among those with glucose
qs_cv <- try(quantile(dm28$gluc_cv, c(1/3,2/3), na.rm=TRUE), silent=TRUE)
qs_mn <- try(quantile(dm28$gluc_mean, c(1/3,2/3), na.rm=TRUE), silent=TRUE)
strat_hr <- list()
if (!inherits(qs_cv,"try-error")) {
  dm28$cv_ter <- cut(dm28$gluc_cv, breaks=c(-Inf, qs_cv, Inf), labels=c("CVt1_low","CVt2","CVt3_high"))
  for (lv in levels(dm28$cv_ter)) {
    strat_hr[[length(strat_hr)+1]] <- fit_tvc_cox(dm28[!is.na(dm28$cv_ter) & dm28$cv_ter==lv,],
      "time_28","status_28","death_28d", paste("DM", lv))
  }
}
if (!inherits(qs_mn,"try-error")) {
  dm28$mn_ter <- cut(dm28$gluc_mean, breaks=c(-Inf, qs_mn, Inf), labels=c("meant1_low","meant2","meant3_high"))
  for (lv in levels(dm28$mn_ter)) {
    strat_hr[[length(strat_hr)+1]] <- fit_tvc_cox(dm28[!is.na(dm28$mn_ter) & dm28$mn_ter==lv,],
      "time_28","status_28","death_28d", paste("DM", lv))
  }
}
for (hv in c(0,1)) {
  strat_hr[[length(strat_hr)+1]] <- fit_tvc_cox(dm28[dm28$hypo70==hv,],
    "time_28","status_28","death_28d", paste0("DM hypo70=", hv))
}
# SHR subset exploratory
shr_dat <- dm28[!is.na(dm28$shr), ]
if (nrow(shr_dat)>50) {
  strat_hr[[length(strat_hr)+1]] <- fit_tvc_cox(shr_dat, "time_28","status_28","death_28d","DM SHR subset (has HbA1c)")
}
strat_bind <- do.call(rbind, strat_hr[!sapply(strat_hr,is.null)])
if (!is.null(strat_bind) && nrow(strat_bind)) write_both(strat_bind, "10_hr_dm_glucose_strata.csv", c("n","n_event"))

# Secondary: stroke logistic; discharge AF logistic; ICU LOS Fine-Gray
log("secondary outcomes")
sec <- list()
# stroke
dsec <- dat[in3p, ]
dsec$exp <- factor(dsec$af_group, levels=c("no_AF","NOAF","preexisting"))
for (dm_lab in c("DM","nonDM","all")) {
  s <- dsec
  if (dm_lab=="DM") s <- s[s$has_dm==1,]
  if (dm_lab=="nonDM") s <- s[s$has_dm==0,]
  s2 <- s[stats::complete.cases(s[, c("has_stroke","exp", covars)]), ]
  if (nrow(s2)<50 || sum(s2$has_stroke)<10) next
  fit <- try(glm(as.formula(paste("has_stroke ~ exp +", paste(covars, collapse="+"))), data=s2, family=binomial()), silent=TRUE)
  if (inherits(fit,"try-error")) next
  sm <- summary(fit)$coefficients; ci <- confint.default(fit)
  for (tm in grep("^exp", rownames(sm), value=TRUE)) {
    sec[[length(sec)+1]] <- data.frame(model="logistic_cc", outcome="stroke_ICD", subset=dm_lab, term=sub("^exp","",tm),
      hr=exp(sm[tm,1]), lcl=exp(ci[tm,1]), ucl=exp(ci[tm,2]), p=sm[tm,4], n=nrow(s2), n_event=sum(s2$has_stroke), stringsAsFactors=FALSE)
  }
}
# discharge AF among ICU survivors with last rhythm
dsurv <- dat[in3p & dat$died_in_icu==0 & !is.na(dat$last_is_af), ]
dsurv$exp <- factor(dsurv$af_group, levels=c("no_AF","NOAF","preexisting"))
dsurv$dc_af <- as.integer(dsurv$last_is_af==1)
for (dm_lab in c("DM","all")) {
  s <- if (dm_lab=="DM") dsurv[dsurv$has_dm==1,] else dsurv
  s2 <- s[stats::complete.cases(s[, c("dc_af","exp", covars)]), ]
  if (nrow(s2)<50) next
  # no_AF has essentially no discharge AF by construction -> complete separation
  tab <- table(s2$exp, s2$dc_af)
  if (any(tab==0) || min(rowSums(tab))<5) {
    log("skip discharge_AF logistic (separation or sparse):", dm_lab)
    next
  }
  fit <- try(glm(as.formula(paste("dc_af ~ exp +", paste(covars, collapse="+"))), data=s2, family=binomial()), silent=TRUE)
  if (inherits(fit,"try-error")) next
  sm <- summary(fit)$coefficients; ci <- confint.default(fit)
  for (tm in grep("^exp", rownames(sm), value=TRUE)) {
    orv <- exp(sm[tm,1]); lcl <- exp(ci[tm,1]); ucl <- exp(ci[tm,2])
    if (!is.finite(orv) || !is.finite(lcl) || !is.finite(ucl) || orv>100 || ucl>1e6) {
      log("skip discharge_AF term (non-finite/separated):", dm_lab, tm)
      next
    }
    sec[[length(sec)+1]] <- data.frame(model="logistic_cc", outcome="discharge_AF", subset=dm_lab, term=sub("^exp","",tm),
      hr=orv, lcl=lcl, ucl=ucl, p=sm[tm,4], n=nrow(s2), n_event=sum(s2$dc_af), stringsAsFactors=FALSE)
  }
}
sec_bind <- do.call(rbind, sec)
if (!is.null(sec_bind) && nrow(sec_bind)) write_both(sec_bind, "11_secondary_logistic.csv", c("n","n_event"))

# Fine-Gray: time to ICU discharge alive, death in ICU competing
if (requireNamespace("cmprsk", quietly=TRUE)) {
  fg_rows <- list()
  for (dm_lab in c("DM","all")) {
    s <- dat[in3p, ]
    if (dm_lab=="DM") s <- s[s$has_dm==1,]
    s$exp <- factor(s$af_group, levels=c("no_AF","NOAF","preexisting"))
    s <- s[stats::complete.cases(s[, c("icu_los_days","died_in_icu","exp")]), ]
    # ftime = icu los; fstatus 1=discharged alive, 2=died in ICU
    s$fstatus <- ifelse(s$died_in_icu==1, 2, 1)
    s$ftime <- pmax(s$icu_los_days, 1/24)
    # unadjusted Fine-Gray for discharge (type 1) by exp; covariates optional via crr but slower
    # crr needs numeric covariates; do unadjusted + age/sofa adjusted if possible
    X <- model.matrix(~ exp + anchor_age + sofa, data=s)[, -1, drop=FALSE]
    ok <- complete.cases(X)
    fit <- try(cmprsk::crr(s$ftime[ok], s$fstatus[ok], X[ok,], failcode=1, cencode=0), silent=TRUE)
    # wait: cencode should be a code not present if everyone has 1 or 2. There is no censor if all leave ICU (dead or alive).
    # All have an event (discharge or death). Use cencode=0 and no zeros in fstatus.
    if (inherits(fit,"try-error")) { log("crr failed", dm_lab, as.character(fit)); next }
    cf <- fit$coef; se <- sqrt(diag(fit$var))
    nm <- names(cf)
    for (i in seq_along(cf)) {
      if (!grepl("^exp", nm[i])) next
      fg_rows[[length(fg_rows)+1]] <- data.frame(model="finegray_crr", outcome="icu_discharge_alive", subset=dm_lab,
        term=sub("^exp","",nm[i]), hr=exp(cf[i]), lcl=exp(cf[i]-1.96*se[i]), ucl=exp(cf[i]+1.96*se[i]),
        p=2*pnorm(-abs(cf[i]/se[i])), n=sum(ok), n_event=sum(s$fstatus[ok]==1), stringsAsFactors=FALSE)
    }
  }
  fg_bind <- do.call(rbind, fg_rows)
  if (!is.null(fg_bind) && nrow(fg_bind)) write_both(fg_bind, "11_finegray_icu_discharge.csv", c("n","n_event"))
}

# ---------------- eICU crude validation ----------------
log("eICU crude validation")
eicu_ok <- try(dbGetQuery(con, "SELECT 1 FROM wh.eicu_parquet.patient LIMIT 1"), silent=TRUE)
if (!inherits(eicu_ok, "try-error")) {
  dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE eicu_first AS
  WITH first AS (
    SELECT patientunitstayid, patienthealthsystemstayid, uniquepid,
           age, gender, hospitaldischargestatus, unitdischargestatus,
           hospitaladmitoffset, unitdischargeoffset, hospitaldischargeoffset
    FROM wh.eicu_parquet.patient
    QUALIFY ROW_NUMBER() OVER (PARTITION BY uniquepid ORDER BY patientunitstayid, hospitaladmitoffset) = 1
  )
  SELECT * FROM first
  WHERE TRY_CAST(age AS INTEGER) >= 18 OR age='> 89'
  ")
  # AF / DM from diagnosis and pastHistory (same spirit as feasibility)
  dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE eicu_dx AS
  SELECT patientunitstayid,
    MIN(diagnosisoffset) FILTER (WHERE regexp_matches(lower(coalesce(diagnosisstring,'')||' '||coalesce(icd9code,'')),
        'atrial fib|atrial flutter|a-fib|\\\\bi48|427.31|42731')) AS af_offset,
    MAX(CASE WHEN regexp_matches(lower(coalesce(diagnosisstring,'')||' '||coalesce(icd9code,'')),
        'atrial fib|atrial flutter|a-fib|\\\\bi48|427.31|42731') THEN 1 ELSE 0 END) AS has_af_dx,
    MAX(CASE WHEN regexp_matches(lower(coalesce(diagnosisstring,'')||' '||coalesce(icd9code,'')),
        'diabetes|\\\\be1[0-4]|\\\\b250') THEN 1 ELSE 0 END) AS has_dm_dx
  FROM wh.eicu_parquet.diagnosis
  GROUP BY 1
  ")
  dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE eicu_ph AS
  SELECT patientunitstayid,
    MAX(CASE WHEN regexp_matches(lower(coalesce(pasthistoryvalue,'')||' '||coalesce(pasthistorypath,'')),
        'atrial fib|atrial flutter|a-fib') THEN 1 ELSE 0 END) AS ph_af,
    MAX(CASE WHEN regexp_matches(lower(coalesce(pasthistoryvalue,'')||' '||coalesce(pasthistorypath,'')),
        'diabetes') THEN 1 ELSE 0 END) AS ph_dm
  FROM wh.eicu_parquet.pastHistory
  GROUP BY 1
  ")
  dbExecute(con, "
  CREATE OR REPLACE TEMP TABLE eicu_class AS
  SELECT
    f.*,
    CASE WHEN coalesce(p.ph_dm,0)=1 OR coalesce(d.has_dm_dx,0)=1 THEN 1 ELSE 0 END AS has_dm,
    CASE WHEN coalesce(p.ph_af,0)=1 OR (coalesce(d.has_af_dx,0)=1 AND d.af_offset<=120) THEN 'preexisting'
         WHEN coalesce(p.ph_af,0)=0 AND coalesce(d.has_af_dx,0)=1 AND d.af_offset>120 THEN 'possible_NOAF'
         WHEN coalesce(p.ph_af,0)=0 AND coalesce(d.has_af_dx,0)=0 THEN 'no_AF'
         ELSE 'indeterminate' END AS eicu_group,
    CASE WHEN lower(coalesce(f.hospitaldischargestatus,'')) IN ('expired','dead') THEN 1 ELSE 0 END AS hosp_death
  FROM eicu_first f
  LEFT JOIN eicu_dx d USING (patientunitstayid)
  LEFT JOIN eicu_ph p USING (patientunitstayid)
  ")
  eicu_tab <- dbGetQuery(con, "
  SELECT has_dm, eicu_group, COUNT(*) AS n,
         COUNT(*) FILTER (WHERE hosp_death=1) AS n_hosp_death
  FROM eicu_class
  GROUP BY 1,2
  ORDER BY 1 DESC, 2
  ")
  eicu_tab$pct_hosp_death <- ifelse(eicu_tab$n>0, round(100*eicu_tab$n_hosp_death/eicu_tab$n, 1), NA)
  write_both(eicu_tab, "12_eicu_crude_counts.csv", c("n","n_hosp_death"))
  # crude logistic NOAF vs no_AF, adjust age/sex if possible
  eicu_df <- dbGetQuery(con, "
  SELECT has_dm, eicu_group, hosp_death,
         CASE WHEN age='> 89' THEN 90 ELSE TRY_CAST(age AS INTEGER) END AS age_n,
         gender
  FROM eicu_class
  WHERE eicu_group IN ('possible_NOAF','no_AF','preexisting')
  ")
  eicu_df$exp <- factor(eicu_df$eicu_group, levels=c("no_AF","possible_NOAF","preexisting"))
  eicu_or <- list()
  for (dm_lab in c("DM","nonDM","all")) {
    s <- eicu_df
    if (dm_lab=="DM") s <- s[s$has_dm==1,]
    if (dm_lab=="nonDM") s <- s[s$has_dm==0,]
    s <- s[stats::complete.cases(s[, c("hosp_death","exp","age_n")]), ]
    if (nrow(s)<50 || sum(s$hosp_death)<10) next
    fit <- try(glm(hosp_death ~ exp + age_n + gender, data=s, family=binomial()), silent=TRUE)
    if (inherits(fit,"try-error")) next
    sm <- summary(fit)$coefficients; ci <- confint.default(fit)
    for (tm in grep("^exp", rownames(sm), value=TRUE)) {
      eicu_or[[length(eicu_or)+1]] <- data.frame(subset=dm_lab, term=sub("^exp","",tm),
        or=exp(sm[tm,1]), lcl=exp(ci[tm,1]), ucl=exp(ci[tm,2]), p=sm[tm,4],
        n=nrow(s), n_event=sum(s$hosp_death), stringsAsFactors=FALSE)
    }
  }
  eicu_or_b <- do.call(rbind, eicu_or)
  if (!is.null(eicu_or_b) && nrow(eicu_or_b)) write_both(eicu_or_b, "12_eicu_crude_or.csv", c("n","n_event"))
  log("eICU done")
} else {
  log("eICU views not available; skip")
}

# ---------------- Forest plots ----------------
log("forest plots")
plot_forest <- function(hrdf, title, fname, term_filter="NOAF") {
  d <- hrdf[!is.na(hrdf$hr) & grepl(term_filter, hrdf$term), ]
  d <- d[is.finite(d$hr) & is.finite(d$lcl) & is.finite(d$ucl) & d$hr>0 & d$lcl>0 & d$ucl<1e6, ]
  if (!nrow(d)) return(invisible(NULL))
  d$lab <- paste0(d$subset, " | ", d$outcome, " | ", d$model)
  d <- d[order(d$subset, d$outcome, d$model), ]
  d$lab <- factor(d$lab, levels=rev(unique(d$lab)))
  p <- ggplot(d, aes(x=hr, y=lab)) +
    geom_vline(xintercept=1, linetype=2, colour="grey40") +
    geom_errorbar(aes(xmin=lcl, xmax=ucl), orientation="y", width=0.22, linewidth=0.4) +
    geom_point(size=2, colour=okabe[["NOAF"]]) +
    scale_x_log10() +
    labs(x="Hazard ratio (95% CI) vs no AF", y=NULL, title=title,
         caption="Time-varying Cox unless otherwise labelled. SAP v1.0.") +
    theme_journal(10) + theme(panel.grid.major.y=element_blank())
  save_fig(p, fname, width=8.8, height=max(4.5, 0.28*nrow(d)+2))
}
if (exists("hr_all") && !is.null(hr_all) && nrow(hr_all)) {
  plot_forest(hr_all[hr_all$model %in% c("tvc_complete_case","tvc_mice_m10","logistic_complete_case"), ],
              "Main models: NOAF vs no AF", "fig3_forest_main")
}
if (exists("sens_bind") && !is.null(sens_bind) && nrow(sens_bind)) {
  plot_forest(sens_bind, "Sensitivity analyses (28-day death): NOAF vs no AF", "fig3_forest_sensitivity")
}
if (exists("strat_bind") && !is.null(strat_bind) && nrow(strat_bind)) {
  plot_forest(strat_bind, "Diabetes glucose strata (28-day death): NOAF vs no AF", "fig4_forest_glucose_strata")
}

# ---------------- Session info ----------------
sess <- data.frame(
  item=c("date","timezone","R","duckdb","survival","mice","seed","SAP"),
  value=c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Pacific/Apia (UTC+13)",
          R.version.string, as.character(packageVersion("duckdb")),
          as.character(packageVersion("survival")),
          if (requireNamespace("mice", quietly=TRUE)) as.character(packageVersion("mice")) else "not used",
          as.character(SEED), "docs/SAP_v1.md"),
  stringsAsFactors=FALSE
)
write.csv(sess, file.path(out_dir, "00_session_info.csv"), row.names=FALSE, fileEncoding="UTF-8")
log("DONE main analysis")
