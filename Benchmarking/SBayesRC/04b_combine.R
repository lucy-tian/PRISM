#!/usr/bin/env Rscript

# threads from Slurm
thr <- Sys.getenv("SLURM_CPUS_PER_TASK", "4")
Sys.setenv(OMP_NUM_THREADS = thr, OPENBLAS_NUM_THREADS = thr, MKL_NUM_THREADS = thr)

suppressPackageStartupMessages({
  library(SBayesRC)
})

# ---------- args ----------
args  <- commandArgs(trailingOnly = TRUE)
trait <- args[[1]]

# ---------- PRS files (3 cols: FID IID PRS, no header) ----------
in_dir  <- "prs"
eur_prs <- file.path(in_dir, paste0("WB_",  trait, ".score.txt"))  # WB/EUR PRS on AFR-val+test
afr_prs <- file.path(in_dir, paste0("Afr_", trait, ".score.txt"))  # AFR PRS on same samples

# ---------- tuning IDs (2 cols: FID IID, no header) ----------
tune_id <- "Afr_val.keep"

# ---------- phenotype for tuning IDs only (3 cols: FID IID PHENO, no header) ----------
pheno_dir  <- "pheno"
pheno_file <- file.path(pheno_dir, paste0("test_val_", trait, ".pheno.txt"))

# ---------- output ----------
out_dir   <- "prs"
outPrefix <- file.path(out_dir, paste0("tuned_afr_wb_", trait))

# (optional) quick sanity checks
stopifnot(file.exists(eur_prs), file.exists(afr_prs), file.exists(tune_id), file.exists(pheno_file))

# ---------- run ----------
SBayesRC::sbrcMulti(
  prs1     = eur_prs,
  prs2     = afr_prs,
  outPrefix= outPrefix,
  tuneid   = tune_id,
  pheno    = pheno_file
)
