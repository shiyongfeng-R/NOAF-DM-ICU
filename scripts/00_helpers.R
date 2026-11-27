# Shared helpers. Sourced by analysis scripts. Do not print individual-level data.
suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(survival); library(ggplot2)
})

# Path configuration: no personal home or machine-specific paths are hard-coded.
# Required environment variables (see README.md):
#   NOAF_ROOT          Warehouse root that contains catalog/coredata.duckdb
#   NOAF_MIMIC_DUCKDB  Path to mimic.duckdb (attached READ_ONLY)
# Optional:
#   NOAF_PROJECT       Directory for output/, output_addon/, logs/
#                      (default: file.path(NOAF_ROOT, "projects/af_dm_noaf_main"))
#   NOAF_REPO          Root of this cloned repository (auto-detected if unset)
#   NOAF_EICU_DUCKDB   Path to eicu.duckdb if a local file is to be attached READ_ONLY
.noaf_need <- function(var) {
  v <- Sys.getenv(var)
  if (!nzchar(v)) stop("Environment variable ", var, " is not set. See README.md.", call. = FALSE)
  path.expand(v)
}
.noaf_detect_repo <- function() {
  v <- Sys.getenv("NOAF_REPO")
  if (nzchar(v) && file.exists(file.path(path.expand(v), "scripts", "00_helpers.R"))) {
    return(normalizePath(path.expand(v)))
  }
  if (file.exists(file.path("scripts", "00_helpers.R"))) return(normalizePath("."))
  if (file.exists("00_helpers.R")) return(normalizePath(".."))
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(f)) {
    cand <- dirname(dirname(normalizePath(f[1])))
    if (file.exists(file.path(cand, "scripts", "00_helpers.R"))) return(cand)
    cand2 <- dirname(cand)
    if (file.exists(file.path(cand2, "scripts", "00_helpers.R"))) return(cand2)
  }
  stop("Could not locate this repository. Set NOAF_REPO to the cloned repo root.", call. = FALSE)
}
repo <- .noaf_detect_repo()
root <- Sys.getenv("NOAF_ROOT")
proj <- Sys.getenv("NOAF_PROJECT")
if (nzchar(root)) root <- path.expand(root)
if (nzchar(proj)) {
  proj <- path.expand(proj)
} else {
  if (!nzchar(root)) stop("Set NOAF_ROOT or NOAF_PROJECT (see README.md).", call. = FALSE)
  proj <- file.path(root, "projects/af_dm_noaf_main")
}
mim_duckdb <- Sys.getenv("NOAF_MIMIC_DUCKDB")
if (nzchar(mim_duckdb)) mim_duckdb <- path.expand(mim_duckdb)
eicu_duckdb <- Sys.getenv("NOAF_EICU_DUCKDB")
if (nzchar(eicu_duckdb)) eicu_duckdb <- path.expand(eicu_duckdb)
out_dir <- file.path(proj, "output")
int_dir <- file.path(out_dir, "INTERNAL_ONLY")
fig_dir <- file.path(out_dir, "figures")
log_dir <- file.path(proj, "logs")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(int_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

SEED <- 20261120L
logf <- file.path(log_dir, "01_run_main.log")
log <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse=" "))
  cat(msg, "\n", sep=""); cat(msg, "\n", file=logf, append=TRUE)
}

exists_ok <- function(path) file.exists(path) && isTRUE(file.info(path)$size > 20)

write_internal <- function(df, name) {
  write.csv(df, file.path(int_dir, name), row.names=FALSE, fileEncoding="UTF-8", na="")
}
write_public <- function(df, name, count_cols=NULL) {
  pub <- df
  if (is.null(count_cols)) {
    count_cols <- names(pub)[grepl("^(n_|n$)|_n$", names(pub)) | names(pub) %in% c("n","N","Freq")]
  }
  n_small <- function(x) { x <- as.numeric(x); !is.na(x) & x <= 9 }
  small_by_col <- list()
  row_denom_small <- rep(FALSE, nrow(df))
  for (cc in count_cols) {
    if (!cc %in% names(df)) next
    sm <- n_small(df[[cc]])
    small_by_col[[cc]] <- sm
    if (cc %in% c("n","N","n_total","n_event","n_atrisk")) row_denom_small <- row_denom_small | sm
    x <- as.numeric(df[[cc]])
    pub[[cc]] <- ifelse(is.na(x), "", ifelse(sm, "\u22649", format(as.integer(round(x)), scientific=FALSE, trim=TRUE)))
  }
  # Mask a rate only when its matching count or the row denominator is <=9
  rcols <- names(df)[grepl("pct_|rate_|per_100|incidence", names(df), ignore.case=TRUE)]
  for (cc in rcols) {
    mask <- row_denom_small
    cands <- unique(c(sub("^pct_", "n_", cc), sub("_pct$", "", cc), sub("^inc_pct", "n", cc)))
    if (grepl("^pct_miss_", cc)) cands <- c(cands, sub("^pct_miss_", "n_miss_", cc))
    for (cand in cands) if (cand %in% names(small_by_col)) mask <- mask | small_by_col[[cand]]
    pub[[cc]] <- ifelse(mask, "\u22649", as.character(df[[cc]]))
  }
  write.csv(pub, file.path(out_dir, name), row.names=FALSE, fileEncoding="UTF-8", na="")
}
write_both <- function(df, name, count_cols=NULL) {
  write_internal(df, name); write_public(df, name, count_cols)
}

# VanderWeele & Ding 2017 E-value for an HR (or OR treated approximately)
evalue_hr <- function(hr) {
  if (is.null(hr) || length(hr)==0 || is.na(hr) || hr <= 0) return(NA_real_)
  h <- if (hr < 1) 1 / hr else hr
  as.numeric(h + sqrt(h * (h - 1)))
}

okabe <- c(no_AF="#0072B2", NOAF="#D55E00", preexisting="#009E73", grayzone="#999999", indeterminate="#CC79A7")

theme_journal <- function(base=11) {
  theme_bw(base_size=base, base_family="sans") +
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),
      strip.background=element_rect(fill="white", colour="black"),
      legend.position="bottom",
      plot.title=element_text(face="bold", size=base),
      axis.title=element_text(size=base)
    )
}

save_fig <- function(p, name, width=7, height=5.5) {
  pdfp <- file.path(fig_dir, paste0(name, ".pdf"))
  pngp <- file.path(fig_dir, paste0(name, ".png"))
  # cairo needs X11 on this Mac and only warns; require a real file then fall back
  pdf_ok <- FALSE
  try({
    ggsave(pdfp, p, width=width, height=height, device=cairo_pdf)
    pdf_ok <- isTRUE(file.exists(pdfp) && file.info(pdfp)$size > 20)
  }, silent=TRUE)
  if (!pdf_ok) {
    if (file.exists(pdfp)) unlink(pdfp)
    ggsave(pdfp, p, width=width, height=height, device=grDevices::pdf, useDingbats=FALSE)
  }
  ggsave(pngp, p, width=width, height=height, dpi=300, device="png")
  log("figure", name, "->", pdfp, pngp)
  invisible(list(pdf=pdfp, png=pngp))
}

year_primary_levels <- c("2008 - 2010", "2011 - 2013", "2014 - 2016", "2017 - 2019")
year_all_levels <- c(year_primary_levels, "2020 - 2022")
