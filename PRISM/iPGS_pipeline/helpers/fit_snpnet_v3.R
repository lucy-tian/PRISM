fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

####################################################################
# command line args
####################################################################

save_d           <- args[1]
phenotype_file   <- args[2]
keep_file        <- args[3]  # "None" or path to .keep file
phenotype        <- args[4]
family           <- args[5]
refit_rdata_file <- args[6]  # "None" or path to .RData file
genotype_pfile   <- args[7]
covariates_type  <- args[8]
p_factor_file    <- args[9]  # "None" or path to .rds file
sample_wts_file  <- args[10] # "None" or path to .rds file

####################################################################
# check the input args
####################################################################
check_file_exists_or_none <- function(f) {
  stopifnot(f == "None" | file.exists(f))
}

stopifnot(dir.exists(save_d))
stopifnot(file.exists(phenotype_file))
check_file_exists_or_none(keep_file)
stopifnot(family %in% c("gaussian", "binomial"))
check_file_exists_or_none(refit_rdata_file)
for (ext in c("pgen", "psam", "pvar.zst")) {
  stopifnot(file.exists(sprintf("%s.%s", genotype_pfile, ext)))
}
check_file_exists_or_none(p_factor_file)
#stopifnot(covariates_type %in% c("UKB_18PCs", "pop_10PCs", "hgdpKGP_10PCs"))
check_file_exists_or_none(sample_wts_file)

####################################################################
# paths and functions
####################################################################

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

script_name <- normalizePath(
    sub("--file=", "", fullargs[grep("--file=", fullargs)])
)
script_dir <- dirname(script_name)

source(file.path(dirname(script_dir), "helpers", "functions.R"))
source(file.path(dirname(script_dir), "paths.sh"))

# For the path to a given RData file
# (saved from a snpnet run with validation set),
# load the lambda sequence stored in the RData file and
# return the lambda sequence to the best lambda value
# (based on the validation set metric).
# This is useful for performing a "refit" using
# a combined set of training + validation set.
read_lambda_sequence <- function(rdata_file) {
  new_env <- new.env()
  load(rdata_file, envir = new_env)
  lambda <- with(
    new_env$fit_snpnet, full.lams[1:which.max(metric.val)]
  )
  rm(new_env)
  return(lambda)
}

####################################################################
# packages
####################################################################

library(snpnet)
configs <- list(
    
  plink2.path = "plink2",   # path to plink2 program
  zstdcat.path = "zstdcat"  # path to zstdcat program
)
# check if the provided paths are valid
for (name in names(configs)) {
  tryCatch(system(paste(configs[[name]], "-h"), ignore.stdout = T),
    condition = function(e) cat("Please add", configs[[name]], "to PATH, or modify the path in the configs list.")
  )
}
# devtools::load_all('/net/bmc-lab5/data/kellis/users/tanigawa/repos/rivas-lab/snpnet')
#devtools::load_all(file.path(save_d, "snpnet"))

####################################################################
# file paths and parameters
####################################################################

covariates <- get_covariates(covariates_type)
is_refit <- (refit_rdata_file != "None")

if (is_refit) {
  split.col <- NULL
  lambda <- read_lambda_sequence(refit_rdata_file)
  niter <- 50 # default value
}else{
  split.col <- "split"
  lambda <- NULL
  niter <- 300
}

if (p_factor_file == "None") {
  p_factor <- NULL
} else {
  p_factor <- read_pfactor(p_factor_file, genotype_pfile, length(covariates))
}

if (sample_wts_file == "None") {
  weights <- NULL
} else {
  weights <- readRDS(sample_wts_file)
}

####################################################################

# call snpnet:snpnet()
#weights = weights
fit_snpnet <- snpnet(
  genotype.pfile = genotype_pfile,
  phenotype.file = phenotype_file,
  family = family,
  phenotype = phenotype,
  covariates = covariates,
  split.col = split.col,
  p.factor = p_factor,
  mem = 60000,  # amount of memory (MB)
  nlambda = 300,
  alpha = .99, # Elastic net penalty
  configs = list(
    results.dir = save_d,  # needed when saving intermediate results
    save = TRUE,  # save intermediate results per iteration (default FALSE)
    nCores = 12,  # number of cores available (default 1)
    keep = keep_file,
    num.snps.batch = 5000,
    lambda = lambda,
    niter = niter,  # max number of iterations (default 50)
    nlams.init = 50,
    prevIter = find_prevIter(save_d),
    verbose = FALSE,
    KKT.verbose = FALSE,
    KKT.thresh = 0,
    plink2.path = "plink2",   # path to plink2 program
    zstdcat.path = "zstdcat"  # path to zstdcat program
  )
)

if (is_refit || (find_prevIter(save_d) < niter)) {
  # check if this is a refit run or we had an early termination

  # save the results
  save(fit_snpnet, file = file.path(save_d, paste0("snpnet.RData")))

  # remove intermediate files
  for (sub_d in c("meta", "results")) {
  if (dir.exists(file.path(save_d, sub_d))) {
    system(sprintf("rm -rf %s", file.path(save_d, sub_d)))
  }
  }

  # export metrics
  get_snpnet_fit_metrics(fit_snpnet) %>%
  rename("#idx" = "idx") %>%
  fwrite(
    file.path(save_d, "snpnet.metrics.tsv.gz"),
    sep = "\t", na = "NA", quote = F
  )

  # export BETAs for idx_export
  idx_export <- with(
      fit_snpnet,
      ifelse(all(is.na(metric.val)), length(beta), which.max(metric.val))
  )
  
  get_snpnet_beta(fit_snpnet, idx_export) %>%
  rename("#ID" = "ID") %>%
  fwrite(
    file.path(save_d, "snpnet.BETAs.tsv.gz"),
    sep = "\t", na = "NA", quote = F
  )
}
