#!/usr/bin/env Rscript

# threads from Slurm
thr <- Sys.getenv("SLURM_CPUS_PER_TASK", "4")
Sys.setenv(OMP_NUM_THREADS = thr, OPENBLAS_NUM_THREADS = thr, MKL_NUM_THREADS = thr)

suppressPackageStartupMessages({
  library(SBayesRC)
})

# ---------- args ----------
args  <- commandArgs(trailingOnly = TRUE)
pop   <- args[[1]]
trait <- args[[2]]

# ---------- weights file ----------
in_dir <- "out"
w <- file.path(in_dir, paste0(pop, "_", trait, "_sbrc.txt"))

# ---------- genotype (PLINK1 BED prefix) ----------
geno_dir    <- "genotype"
geno_prefix <- file.path(geno_dir, "Afr_val_test")  # <- no extra ')'
geno_chr    <- ""                               # empty means single prefix (not per-chr)

# ---------- out ----------
out_dir    <- "prs"
out_prefix <- file.path(out_dir, paste0(pop, "_", trait))

# optional sanity checks (comment out if you prefer)
if (!file.exists(w)) stop("Weight file not found: ", w)
if (!all(file.exists(paste0(geno_prefix, c(".bed",".bim",".fam")))))
  stop("Genotype BED/BIM/FAM not found for prefix: ", geno_prefix)

# run
SBayesRC::prs(weight = w, genoPrefix = geno_prefix, genoCHR = geno_chr, out = out_prefix)
