#!/usr/bin/env Rscript

# threads from Slurm
thr <- Sys.getenv("SLURM_CPUS_PER_TASK", "4")
Sys.setenv(OMP_NUM_THREADS=thr, OPENBLAS_NUM_THREADS=thr, MKL_NUM_THREADS=thr)

suppressPackageStartupMessages({
  library(SBayesRC)
})

# ---------- args ----------
args <- commandArgs(trailingOnly = TRUE)
pop   <- args[[1]]
trait <- args[[2]]

# ---------- paths ----------
annot_dir <- "/home/lucytian/lab4/sbayes/annot_baseline2.2.txt"   # <- use the **unzipped** folder (not .zip)
ld_WB     <- "/home/lucytian/lab4/sbayes/ukbEUR_Imputed"      # must contain snp.info
ld_AFR    <- "/home/lucytian/lab4/sbayes/ukbAFR_Imputed"
sum_in    <- "input_sumstat"
out_dir   <- "out"
prs_dir   <- "prs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(prs_dir, showWarnings = FALSE, recursive = TRUE)

# choose LD dir by pop
ld_dir <- if (pop == "WB") ld_WB else ld_AFR

# input/output prefixes
ma_in   <- file.path(sum_in, sprintf("%s_%s.ma", pop, trait))   # e.g. WB_INI30120.ma
out_pre <- file.path(out_dir, sprintf("%s_%s", pop, trait))     # e.g. out/WB_INI30120


# ---------- run ----------
SBayesRC::tidy(
  mafile = ma_in, LDdir = ld_dir,
  output = paste0(out_pre, "_tidy.ma"), log2file = TRUE
)

SBayesRC::impute(
  mafile = paste0(out_pre, "_tidy.ma"), LDdir = ld_dir,
  output = paste0(out_pre, "_imp.ma"), log2file = TRUE
)

SBayesRC::sbayesrc(
  mafile    = paste0(out_pre, "_imp.ma"),
  LDdir     = ld_dir,
  outPrefix = paste0(out_pre, "_sbrc"),
  annot     = annot_dir,
  log2file  = TRUE,
  tuneStep = c(0.995, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3)
)

cat("Done:", pop, trait, "\n")
