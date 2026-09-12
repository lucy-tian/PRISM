fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

save_d <- args[1]
q <- as.numeric(args[2])

if (!dir.exists(save_d)) dir.create(save_d, recursive = TRUE, showWarnings = FALSE)

library(snpnet)
library(tidyverse)
library(data.table)

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

configs <- list(
  plink2.path = "plink2",   # path to plink2 program
  zstdcat.path = "zstdcat"  # path to zstdcat program
)

for (name in names(configs)) {
  tryCatch(system(paste(configs[[name]], "-h"), ignore.stdout = T),
    condition = function(e) cat("Please add", configs[[name]], "to PATH, or modify the path in the configs list.")
  )
}


is_refit <- FALSE
niter <- 300
genotype_pfile <- "synthetic_v1_Afr_chr22only"
phenotype_file <- "final_pheno.pheno"
keep_file <- "Afr_train_val.keep.list"
p_factor_file <- file.path("penalty_files", paste0("set_P_q", as.integer(q * 100), ".rds"))


read_pfactor <- function(p.factor.file, genotype.pfile, covar_length = 0){
    pvar_variants <- fread(
        cmd = sprintf("zstdcat %s", sprintf("%s.pvar.zst", genotype.pfile)),
        select = c("ID", "ALT")
    ) %>%
    mutate(ID_ALT = paste0(ID, "_", ALT)) %>%
    pull(ID_ALT)

    p.factor <- readRDS(p.factor.file)[pvar_variants]

    # normalize p.factor so that the sum is equal to the number of
    # features (genetic variants + covariates)
    p.factor <- (covar_length + length(p.factor)) * p.factor / sum(p.factor)

    return(p.factor)
}

find_prevIter <- function(save_d){
    fs <- Sys.glob(file.path(save_d, "results", "output_iter_*.RData"))

    if(length(fs) == 0){
        return(0)
    }else{
        sapply(fs, function(f){as.integer(str_replace_all(basename(f), "^output_iter_|.RData$", ""))}) %>%
        unname() %>%
        sort() %>%
        last(1) %>%
        return()
    }
}


p_factor <- read_pfactor(p_factor_file, genotype_pfile, 10)

fit_snpnet <- snpnet(
  genotype.pfile = genotype_pfile,
  phenotype.file = phenotype_file,
  family = "gaussian",
  phenotype = "y",
  covariates = paste0("PC", 1:10),
  split.col = "split",
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
    lambda = NULL,
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
}
