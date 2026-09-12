#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# 1. Setup paths and variables
# Usage: Rscript eval_script.R [save_directory]
if (length(args) < 1) {
    stop("Usage: Rscript eval_script.R <save_directory>")
}

save_d     <- args[1]
pgs_method <- "snpnet" 
phenotype  <- "y"       

library(tidyverse)
library(data.table)
library(psychometric)

eval_f <- file.path(save_d, sprintf("%s.eval", pgs_method))
phenotype_file <- "final_pheno.pheno"
sscore_f <- file.path(save_d, "snpnet.sscore.zst") 

# 2. Load and merge data
if (!file.exists(phenotype_file)) { stop(paste("Phenotype file not found:", phenotype_file)) }
if (!file.exists(sscore_f)) { stop(paste("Score file not found:", sscore_f)) }

phe_df <- fread(phenotype_file) %>%
    rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
    rename(population = POP) %>% 
    mutate(across(all_of(c("FID", "IID")), as.character))

# Read .zst file using zstdcat
sscore_df <- fread(cmd = paste("zstdcat", sscore_f)) %>%
    rename_with(~str_replace(., "#", ""), starts_with("#")) %>%
    mutate(across(all_of(c("FID", "IID")), as.character))

# Identify the score column (Folder_SUM or SCORE1_SUM)
sscore_f_col_name <- sprintf("%s_SUM", basename(save_d))
actual_score_col <- intersect(names(sscore_df), c(sscore_f_col_name, "SCORE1_SUM"))[1]

if(is.na(actual_score_col)) {
    actual_score_col <- names(sscore_df)[grep("_SUM$", names(sscore_df))][1]
}

if(is.na(actual_score_col)) {
    stop("Could not find a SUM score column in the .sscore.zst file!")
}

col_pgs <- paste0("PGS_", phenotype)
pc_cols <- paste0("PC", 1:10)

all_df <- phe_df %>%
    inner_join(
        sscore_df %>% dplyr::select(FID, IID, !!sym(actual_score_col)), 
        by = c("FID", "IID")
    ) %>%
    rename(!!col_pgs := !!sym(actual_score_col))

cat(sprintf("Merged data contains %d individuals.\n", nrow(all_df)))

# --- Helper Functions ---

fit_glm <- function(data_df, formula_str, family){
    glm(stats::as.formula(formula_str), family = family, data = data_df)
}

fit_to_df <- function(fit){
    fit_df <- summary(fit)$coeff %>%
        as.data.frame() %>% 
        rownames_to_column("variable")
    colnames(fit_df) <- c("variable", "estimate", "SE", "z_or_t_value", "P")
    return(fit_df)
}

glm_fit_to_R2 <- function(glm_fit) {
    with(summary(glm_fit), 1 - deviance / null.deviance)
}

compose_regression_formula_str <- function(response, predictors, quote_char="`") {
    return(sprintf(
        "%s ~ 1 + %s",
        paste0(quote_char, response, quote_char),
        paste(sapply(predictors, function(term){paste0(quote_char, term, quote_char)}), collapse = " + ")
    ))
}

eval_R2_CI <- function(data, response, predictors, level=.95) {
    formula_str <- compose_regression_formula_str(response, predictors)
    glm_fit <- fit_glm(data, formula_str, "gaussian")
    
    P_val <- glm_fit %>% 
        fit_to_df() %>%
        dplyr::mutate(variable = str_replace_all(variable, "`", "")) %>%
        dplyr::filter(variable %in% predictors) %>%
        dplyr::pull(P) %>% 
        {if(length(.) > 0) min(.) else NA}
    
    rsq <- glm_fit_to_R2(glm_fit)
    
    psychometric::CI.Rsq(rsq, n=nrow(data), k=length(predictors), level=level) %>%
        dplyr::mutate(metric = "R2", response = response, 
                      predictors = paste(predictors, collapse = "+"),
                      P = P_val, n = nrow(data)) %>%
        dplyr::rename("eval"="Rsq", "l_eval"="LCL", "u_eval"="UCL") %>%
        dplyr::select(response, predictors, metric, eval, l_eval, u_eval, P, n)
}

score_eval_wrapper <- function(df, score_cols, phenotype) {
    target_df <- df %>% dplyr::filter(population == "AFR")
    
    splits <- list(
        list(name = "train_val", values = c("train", "val")),
        list(name = "train",     values = "train"),
        list(name = "val",       values = "val"),
        list(name = "test",      values = "test")
    )
    
    lapply(splits, function(s) {
        sub_df <- target_df %>% dplyr::filter(split %in% s$values)
        if(nrow(sub_df) == 0) return(NULL)
        
        eval_R2_CI(sub_df, phenotype, score_cols) %>%
            dplyr::mutate(split = s$name, population = "AFR")
    }) %>% bind_rows()
}

# 3. Execute 3 Scenarios
eval_df <- bind_rows(
    score_eval_wrapper(all_df, pc_cols, phenotype) %>% mutate(model = "covars"),
    score_eval_wrapper(all_df, col_pgs, phenotype) %>% mutate(model = "PGS"),
    score_eval_wrapper(all_df, c(col_pgs, pc_cols), phenotype) %>% mutate(model = "full")
)

# 4. Save results
eval_df %>%
    dplyr::rename("#response" = "response") %>%
    fwrite(sprintf("%s.tsv.gz", eval_f), sep = "\t", na = "NA", quote = FALSE)

cat("Successfully evaluated all 3 models.\n")
cat("Output saved to:", sprintf("%s.tsv.gz", eval_f), "\n")