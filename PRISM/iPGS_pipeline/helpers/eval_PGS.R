fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

####################################################################
# command line args
####################################################################

save_d    <- args[1]
phenotype <- args[2]
family    <- args[3]
covariates_type <- args[4] # we may able to recover this from the save_d
pgs_method <- args[5] # snpnet PRSCSx
phenotype_file <- args[6]
pheno_names_f <- args[7]

# we will read ${save_d}/${pgs_method}.sscore.zst and
# evaluate against ${phenotype} column in the phenotype_file
# using GLM family (either binomial or gaussian)

stopifnot(family %in% c("gaussian", "binomial"))
#stopifnot(covariates_type %in% c("pop_10PCs", "UKB_18PCs", "hgdpKGP_10PCs"))

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
##CHANGE HERE
eval_f       <- file.path(save_d, sprintf("%s.eval_lead", pgs_method))
pgs_vs_phe_f <- file.path(save_d, sprintf("%s.PGS_vs_phe", pgs_method))

####################################################################
# input files
####################################################################

# Optionally, print the path of the phenotype file if a match was found
if (!is.null(phenotype_file)) {
  cat("Phenotype file path:", phenotype_file, "\n")
}



sscore_f <-
    file.path(save_d, sprintf("%s.sscore.zst", pgs_method))
sscore_f_col_name <-
    sprintf("%s_SUM", basename(save_d))
covar_score_f <- file.path(
    data_d,
    "covars",
    phenotype,
    "covars.score.tsv.gz"
)

#    file.path(
#        '/home/lucytian/data/1_Single_Cell_PRS/0_trial',
#        sprintf(
#           "%s_covars_%s",
#           ifelse(
#                covariates_type == "pop_10PCs",
#               "__POPULATION__", "406k"
#            ),
#           covariates_type
#        ),
#        phenotype,
#        "covars.score.tsv.gz"
#    )

####################################################################
# main
####################################################################

covariates <- get_covariates(covariates_type)

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

sscore_df <-
    fread(cmd = paste("zstdcat", sscore_f)) %>%
    rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
    mutate(across(all_of(c("FID", "IID")), as.character))

if (! covariates_type == "pop_10PCs") {
    covar_score_df <-
        covar_score_f %>%
        fread() %>%
        rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
        mutate(across(all_of(c("FID", "IID")), as.character))

} else {
    # read covariate-only model across populations
    covar_score_df <-
        covar_score_f %>%
        str_replace_all("__POPULATION__", "*") %>%
        Sys.glob() %>%
        sort() %>%
        lapply(function(f) {
            f %>%
            fread() %>%
            rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
            mutate(across(all_of(c("FID", "IID")), as.character))
        }) %>%
        bind_rows()
}


col_covars <- paste0("covars_", phenotype)
col_pgs    <- paste0("PGS_", phenotype)
col_full   <- paste0("full_", phenotype)

all_df <-
    phe_df %>%
    inner_join(
        covar_score_df, by = c("FID", "IID")
    ) %>%
    left_join(
        sscore_df %>%
        select(all_of(c("FID", "IID", sscore_f_col_name))) %>%
        rename(!!col_pgs := 3),
        by = c("FID", "IID")
    ) %>%
    mutate(tmp_col_full := rowSums(across(all_of(c(col_covars, col_pgs))))) %>%
    rename(!!col_full := tmp_col_full)


# evaluate R2 or AUROC
if (! covariates_type == "pop_10PCs") {
    eval_df <-
        bind_rows(
            score_eval_wrapper(all_df, col_covars, phenotype, family, "all") %>%
            mutate(model = "covars"),

            score_eval_wrapper(all_df, col_pgs, phenotype, family, "all") %>%
            mutate(model = "PGS_with_lead"),

            score_eval_wrapper(all_df, col_full, phenotype, family, "all") %>%
            mutate(model = "full_with_lead")
        )
} else {
    eval_df <-
        all_df %>%
        pull(population) %>%
        unique() %>%
        sort() %>%
        lapply(function(p) {
            bind_rows(
                score_eval_wrapper(all_df, col_covars, phenotype, family, p) %>%
                mutate(model = "covars"),

                score_eval_wrapper(all_df, col_pgs, phenotype, family, p) %>%
                mutate(model = "PGS_with_lead"),

                score_eval_wrapper(all_df, col_full, phenotype, family, p) %>%
                mutate(model = "full_with_lead")
            )
        }) %>%
        bind_rows()
}

# save R2 or AUROC
eval_df %>%
rename("#response" = "response") %>%
fwrite(sprintf("%s.tsv.gz", eval_f), sep = "\t", na = "NA", quote = F)

# plot R2 or AUROC
plot_eval <-
    eval_df %>%
    plot_eval_CI(family, phenotype, GBE_names = phe_names_dict) +
    theme(legend.position = "bottom")

for (ext in c("png", "pdf")) {ggsave(
    sprintf("%s.%s", eval_f, ext),
    plot_eval,
    width = 8, height = 8
)}

# prepare data frames for the PGS vs phenotype plots
plot_df <-
    all_df %>%
    filter(split == "test") %>%
    drop_na(all_of(c(phenotype, col_pgs))) %>%
    rename("geno_score" := all_of(col_pgs)) %>%
    rename("phe" := all_of(phenotype))

if (covariates_type == "pop_10PCs") {
    # for plotting purpose, we focus on WB
    plot_df <-
        plot_df %>%
        filter(population == "WB")
}

summary_plot_df <-
    plot_df %>%
    pull(population) %>%
    unique() %>%
    c("all") %>%
    lapply(function(pop) {
        if (family == "gaussian" & pop %in% c("WB", "all")) {
            quantile_bins <- c(0, .0005, .01, (1:19) / 20, .99, .9995, 1)
        } else {
            quantile_bins <- ((0:10) / 10)
        }
        plot_df %>%
            filter(if (pop != "all") population == pop else TRUE) %>%
            mutate(
                geno_score_percentile = rank(geno_score) / n()
            ) %>%
            compute_summary_df(
                "geno_score_percentile", "phe",
                bins = quantile_bins,
                family = family
            ) %>%
            mutate(population = pop)
    }) %>%
    bind_rows()

# save a data frame used for the summary plot
summary_plot_df %>%
rename("#l_bin" = "l_bin") %>%
fwrite(
    sprintf("%s.tsv.gz", pgs_vs_phe_f),
    sep = "\t", na = "NA", quote = F
)

# generate plots
if (family == "gaussian") {
    p1 <-
        plot_df %>%
        plot_PRS_vs_phe() +
        theme(legend.position = c(.1, .8)) +
        labs(
            title = get_plot_label(phenotype, phe_names_dict),
            y = get_plot_label(phenotype, phe_names_dict),
            x = sprintf("%s PGS (Z-score)", pgs_method)
        )

    p2 <- summary_plot_df %>%
        filter(population == "all") %>%
        plot_PRS_bin_vs_phe(mean(plot_df$phe)) +
        labs(
            title = get_plot_label(phenotype, phe_names_dict),
            x = sprintf("%s PGS percentile", pgs_method),
            y = get_plot_label(phenotype, phe_names_dict)
        )

} else if (family == "binomial") {
    p1 <- plot_df %>%
        plot_PRS_binomial() +
        labs(
            title = get_plot_label(phenotype, phe_names_dict),
            x = get_plot_label(phenotype, phe_names_dict),
            y = sprintf("%s PGS (Z-score)", pgs_method)
        )

    p2 <- summary_plot_df %>%
        filter(population == "all") %>%
        plot_PRS_bin_vs_OR() +
        labs(
            title = get_plot_label(phenotype, phe_names_dict),
            x = sprintf("%s PGS percentile", pgs_method)
        )

} else{
    stop(sprintf("%s family is not supported!", family))
}

for (ext in c("png", "pdf")) {ggsave(
    sprintf("%s.%s", pgs_vs_phe_f, ext),
    gridExtra::arrangeGrob(p1, p2, ncol = 2),
    width = 12, height = 6
)}
