fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

####################################################################
# command line args
####################################################################

phenotype       <- args[1]
family          <- args[2]
covariates_type <- args[3]
pop <- args[4]
phenotype_file <- args[5]
pheno_names_f <- args[6]

# we use covariates and evaluate their predictive performance
# against ${phenotype} column in the phenotype_file
# using GLM family (either binomial or gaussian)

stopifnot(family %in% c("gaussian", "binomial"))
#stopifnot(pop %in% c("all", "WB", "NBW", "Afr", "SA"))
# population-specific PCs should NOT be used in conjunction with all

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
source(file.path(dirname(script_dir), "helpers", "eval_functions.R"))
source(file.path(dirname(script_dir), "paths.sh"))

####################################################################
# output files
####################################################################

save_d <- file.path(
    data_d,
    "covars",
    phenotype
)

mkdir_p_recursive(save_d, data_d)
covar_beta_f  <- file.path(save_d, "covars.BETAs.tsv.gz")
covar_score_f <- file.path(save_d, "covars.score.tsv.gz")
covar_eval_f  <- file.path(save_d, "covars.eval")

####################################################################
# input files
####################################################################

#pheno_names_f     <- file.path("/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/pheno", "phe.names.tsv.gz")
#phenotype_file    <- file.path("/net/bmc-lab5/data/kellis/group/tanigawa/share/lucytian/20230626", sprintf(
#    "phe%s.tsv.gz",
#    ifelse(pop == "WB", ".WB", "")
#))


pheno_names_f <- file.path("/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/pheno", "blood_count.counts.tsv.gz")
####################################################################
# main
####################################################################

covariates <- get_covariates(covariates_type)
#APOE_snps <- read.csv("/home/lucytian/group/data/gwas_geno/APOE_covars.csv", header = FALSE, stringsAsFactors = FALSE)[, 1]

phe_names_dict <-
    pheno_names_f %>%
    fread(select = c("#GBE_ID", "GBE_short_name")) %>%
    rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
    deframe()


phe_df <-
    phenotype_file %>%
    fread(
        select = c("#FID", "IID", "population", "split", covariates, phenotype)
    ) %>%
    rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
    mutate(across(all_of(c("FID", "IID")), as.character))

# regression formula
covar_formula_str <-
    sprintf(
        "%s ~ 1 + %s",
        phenotype, paste(covariates, collapse = " + ")
    )

col_covars <- paste0("covars_", phenotype)


# fit covariate-only model using training and validation set individuals
# Assuming 'phenotype' is a variable containing the name of your actual phenotype column
phenotype_column_name <- phenotype

# Now using this variable to dynamically select and manipulate the column in phe_df
if (family == "binomial" && all(unique(phe_df[[phenotype_column_name]]) %in% c(1, 2))) {
    phe_df[[phenotype_column_name]] <- phe_df[[phenotype_column_name]] - 1
}


covar_model_betas_df <-
    phe_df %>%
    filter(split %in% c("train", "val")) %>%
    fit_glm(covar_formula_str, family) %>%
    fit_to_df() %>%
    mutate(population = pop)

covar_model_betas_df %>%
rename("#variable" = "variable") %>%
fwrite(covar_beta_f, sep = "\t", na = "NA", quote = F)


# compute covariate-only score
covar_scores_df <-
    phe_df %>%
    FID_IID_to_rownames() %>%
    compute_matrix_product(
        covar_model_betas_df, covariates,
        c("estimate"), "variable"
    ) %>%
    rename(!!col_covars := 1) %>%
    FID_IID_from_rownames()

# save the individual-level covariate-only score
covar_scores_df %>%
rename("#FID" = "FID") %>%
fwrite(covar_score_f, sep = "\t", na = "NA", quote = F)

# evaluate R2 or AUROC
eval_df <-
    phe_df %>%
    left_join(
        covar_scores_df, by = c("FID", "IID")
    ) %>%
    score_eval_wrapper(col_covars, phenotype, family) %>%
    mutate(model = "covars")

# save R2 or AUROC
eval_df %>%
rename("#response" = "response") %>%
fwrite(sprintf("%s.tsv.gz", covar_eval_f), sep = "\t", na = "NA", quote = F)

# plot R2 or AUROC
plot_eval <-
    eval_df %>%
    plot_eval_CI(family, phenotype, GBE_names = phe_names_dict) +
    theme(legend.position = "none")

for (ext in c("png", "pdf")) {
    ggsave(
        sprintf("%s.%s", covar_eval_f, ext),
        plot_eval,
        width = 8, height = 8
    )
}
