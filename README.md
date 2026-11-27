# NOAF-DM-ICU

Analysis code for the manuscript *New-onset versus pre-existing atrial fibrillation and 28-day mortality in critically ill adults with diabetes: a retrospective cohort study with a time-varying exposure and electrocardiogram cross-checking*.

This repository contains **analysis code and the locked statistical analysis plan only**. It does **not** include MIMIC-IV, MIMIC-IV-ECG, eICU-CRD, derived DuckDB files, or any individual-level or intermediate data.

Planned public location (repository not yet created at the time of manuscript preparation): <https://github.com/shiyongfeng-R/NOAF-DM-ICU>.

## Data access

The analyses use:

- MIMIC-IV version 3.1 (<https://doi.org/10.13026/kpb9-mt58>)
- MIMIC-IV-ECG version 1.0 (<https://doi.org/10.13026/4nqg-sb35>)
- eICU-CRD version 2.0 (<https://doi.org/10.13026/C2WM1R>)

MIMIC-IV and eICU-CRD are available on [PhysioNet](https://physionet.org/) to credentialed users after CITI human-subjects training and a signed data use agreement. MIMIC-IV-ECG is openly available under the Open Database License and was linked to MIMIC-IV under the MIMIC-IV data use agreement. Individual-level data cannot be redistributed by the authors. Users must obtain the databases themselves and must not commit data, credentials, or DuckDB files to this repository.

The scripts attach a local catalog (`catalog/coredata.duckdb` under `NOAF_ROOT`) and `mimic.duckdb` as **READ_ONLY**. They never write to those files.

## How to run

1. Install R 4.5.1 (or a compatible 4.x release) with the packages listed under Session info. The original analysis did not install new packages beyond those already on the authors' machine.
2. Build or point to a DuckDB catalog and a `mimic.duckdb` file derived from your PhysioNet download (READ_ONLY attach).
3. In the shell, set the path variables (do not hard-code personal home paths in the scripts):

```bash
export NOAF_ROOT="/path/to/warehouse"                 # directory that contains catalog/coredata.duckdb
export NOAF_MIMIC_DUCKDB="/path/to/mimic.duckdb"      # READ_ONLY
export NOAF_PROJECT="/path/to/writable/project/dir"   # output/, output_addon/, logs/ are written here
# optional:
# export NOAF_REPO="/path/to/this/clone"
# export NOAF_EICU_DUCKDB="/path/to/eicu.duckdb"
```

`NOAF_PROJECT` defaults to `$NOAF_ROOT/projects/af_dm_noaf_main` if unset. Use a directory **outside** the git clone, or a local folder ignored by `.gitignore`, so that outputs are not committed.

4. From the repository root:

```bash
Rscript --vanilla scripts/01_run_main.R
Rscript --vanilla scripts/B_addon/01_run_addon.R
Rscript --vanilla scripts/B_addon/02_followup_20261122.R
Rscript --vanilla scripts/B_addon/03_landmark_sap_rerun.R
Rscript --vanilla scripts/B_addon/04_fix_addon_output_20261127.R
Rscript --vanilla scripts/B_addon/05_revision_C01_C04.R
Rscript --vanilla scripts/B_addon/06_revision_C04b.R
```

The add-on scripts replay the locked data-build section of `01_run_main.R` in memory and write only aggregate tables (counts of 1–9 are masked as ≤9 in public CSV files). They do not overwrite the locked official `output/` or `results_final/` trees. Scripts do not print individual-level rows.

Random seed: `20261120`.

## Script order

| Order | Script | Role |
|---|---|---|
| 0 | `scripts/00_helpers.R` | Shared paths, logging, public/internal CSV writers, E-value helper |
| 1 | `scripts/01_run_main.R` | Locked SAP v1.0 main analysis (cohort, exposure, imputation, primary Cox models, figures) |
| 2 | `scripts/B_addon/01_run_addon.R` | First add-on battery (competing risks, interaction, landmarks, ECG sub-cohort, etc.) |
| 3 | `scripts/B_addon/02_followup_20261122.R` | Follow-up fixes and additional tables |
| 4 | `scripts/B_addon/03_landmark_sap_rerun.R` | Landmark re-run to the SAP specification |
| 5 | `scripts/B_addon/04_fix_addon_output_20261127.R` | Annotation and output-hygiene fixes (does not change locked estimates) |
| 6 | `scripts/B_addon/05_revision_C01_C04.R` | Revision analyses C01–C04 |
| 7 | `scripts/B_addon/06_revision_C04b.R` | Outcome-conditional (Lash) probabilistic bias analysis, superseding C04 as a correction |

The locked statistical analysis plan is `SAP_v1.md` (v1.0, locked 2026-11-20). In this copy, one home-relative path to a local feasibility report was rephrased so that a personal directory is not published; the scientific content is otherwise unchanged.

## Session info

Recorded on the authors' machine (Pacific/Apia, UTC+13). Seed 20261120 throughout.

- Locked main analysis (`00_session_info`): 2026-11-21; R 4.5.1 (2025-06-13); duckdb 1.5.6; survival 3.8.3; mice 3.18.0.
- Add-on battery (`B00_session_info`, 2026-11-22): the same, plus cmprsk 2.2.12, tidycmprsk 1.1.0, ggplot2 4.0.0, strong 1.0.9; episensr was not installed (the bias analysis is implemented in base R and survival).
- Revision C01–C04 (`C00_revision_session_info`, 2026-11-27): R 4.5.1; survival 3.8.3; mice 3.18.0; duckdb 1.5.6.
- Revision C04b (`C04b_session_info`, 2026-11-27): R 4.5.1; survival 3.8.3; duckdb 1.5.6.

The SAP also lists packages available on the machine but not necessarily loaded in every script: survminer 0.5.1, rms 8.0.0, tableone 0.13.2, forestplot 3.1.7.

## Licence

This code is released under the MIT Licence (see `LICENSE`). The clinical databases remain under their PhysioNet terms and are not included.

## Citation

Please cite the published article when available, and the PhysioNet databases as required by their data use agreements.
