###############################################
# file IO
###############################################

#' mkdir -p
#'
#' @example mkdir_p_recursive(glm_save_d, data_d)
mkdir_p_recursive <- function(dir, root_dir){
    if(file.exists(root_dir)){
        l <- list()
        mkdir_p_dir <- dir
        while(TRUE){
            if(mkdir_p_dir == root_dir){
                break
            }
            l <- c(l, mkdir_p_dir)
            mkdir_p_dir <- dirname(mkdir_p_dir)
        }
        for(d in rev(l)){
            dir.create(d, showWarnings = FALSE)
        }
    }
}


###############################################
# covariates
###############################################

#' Given a covariate type, return the list of covariates
#'
#' @example get_covariates("pop_10PCs")
#' @example get_covariates("UKB_18PCs")
#' @example get_covariates("hgdpKGP_10PCs")


get_covariates <- function(covariates_string) {
    # split on comma, trim whitespace
    covariates <- unlist(strsplit(covariates_string, ","))
    covariates <- trimws(covariates)
    return(covariates)
}

#get_covariates <- function(covariates_type) {
#
#    stopifnot(covariates_type %in% c("pop_10PCs", "UKB_18PCs", "hgdpKGP_10PCs"))
#
#    common_covariates <- c(
#        "age", "sex", "Townsend", "age_age", "age_sex"
#    )
#
#    if (covariates_type == "pop_10PCs") {
#        covariates <- c(
#            common_covariates, paste0("PC", 1:10)
#        )
#
#    } else if (covariates_type == "UKB_18PCs") {
#        covariates <- c(
#            common_covariates, paste0("UKB_PC", 1:18)
#        )
#
#    } else if (covariates_type == "hgdpKGP_10PCs") {
#        covariates <- c(
#            common_covariates, paste0("hgdpKGP_PC", 1:10)
#        )
#    }
#
#    return(covariates)
#}


###############################################
# snpnet
###############################################

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


get_snpnet_fit_metrics <- function(fit_snpnet){
    with(
        fit_snpnet,
        data.frame(
            full.lams = full.lams,
            metric.train = metric.train,
            metric.val = metric.val
        ) %>%
        drop_na(metric.train) %>%
        mutate(
            idx = 1:n(),
            rank.metric.val = rank(-metric.val),
            nzero = sapply(beta, function(b){enframe(b) %>% filter(value != 0) %>% nrow})
        ) %>%
        select(all_of(c(
            "idx", "full.lams", "nzero",
            "metric.train", "metric.val", "rank.metric.val"
        )))
    )
}


get_snpnet_beta <- function(fit_snpnet, lambda_idx) {
    covariates <- fit_snpnet$configs$covariates

    with(
        fit_snpnet,
        beta[[lambda_idx]] %>% enframe(value = "BETA") %>% filter(BETA != 0)
    ) %>%
    separate(
        name, c("ID", "A1"),
        sep = "_", remove = F, extra = "drop", fill = "left"
    ) %>%
    separate(
        ID, c("CHROM", "POS", "REF", "ALT"),
        sep = ":", remove = F, extra = "drop", fill = "left"
    ) %>%
    mutate(
        ID = if_else(name %in% covariates, name, ID),
        A1 = if_else(name %in% covariates, "", A1),
        CHROM = factor(
            if_else(name %in% covariates, "", CHROM),
            levels = c("", 1:22, "X", "Y", "XY", "MT")
        ),
        POS = if_else(name %in% covariates, "", POS),
        REF = if_else(name %in% covariates, "", REF),
        ALT = if_else(name %in% covariates, "", ALT)
    ) %>%
    arrange(CHROM, POS, REF, ALT) %>%
    select(all_of(c("ID", "A1", "BETA")))
}


##################################
#' A wrapper function for glm()
#'
#' Given a data frame, a regression formula in string, and a model family,
#' call glm function.
#'
#' @param data_df A data frame containing the dataset for glm analysis
#' @param formula_str A string representing a regression formula
#' @param family The GLM family
#' @return a glm fit object
#' @examples
#' \dontrun{
#' fit_glm(data, 'response ~ 1 + x + y', "binomial")
#' }
#'
#' @export
fit_glm <- function(data_df, formula_str, family){
    glm(
        stats::as.formula(formula_str),
        family = family,
        data = data_df
    )
}


#' convert the fit object to a data frame
#'
#' @param fit a fit object from glm
#'
#' @export
fit_to_df <- function(fit){
    fit_df <- summary(fit)$coeff %>%
    as.data.frame() %>% rownames_to_column("variable")
    colnames(fit_df) <- c("variable", "estimate", "SE", "z_or_t_value", "P")
    fit_df
}


#' Assign rownames of phenotype data frame based on FID and IID
#'
#' Given a data frame of the phenotype data, concatenate FID and IID
#' columns with a separater character (sep) and use it as a rowname
#'
#' @param pheno_df A phenotype data frame with columns named 'FID' and 'IID'
#' @param sep A separater to concatenate FID and IID
#' @return an updated phenotype data drame with rownames
#' @examples
#' \dontrun{
#' FID_IID_to_rownames(pheno_df)
#' }
#'
#' @export
FID_IID_to_rownames <- function(pheno_df, sep ="_"){
    pheno_df %>%
    mutate(FID_IID = paste(FID, IID, sep = sep)) %>%
    select(-FID, -IID) %>%
    column_to_rownames("FID_IID")
}


#' Recover FID and IID from rownames of the data frame
#'
#' Given a data frame of the phenotype data, recover FID and IID from
#' rownames
#'
#' @param pheno_df A phenotype data frame with rownames
#' @param sep A separater used to concatenate FID and IID
#' @return an updated phenotype data drame with FID and IID columns
#'
#' @examples
#' \dontrun{
#' FID_IID_from_rownames(pheno_df)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#'
#' @export
FID_IID_from_rownames <- function(pheno_df, sep = "_"){
    pheno_df %>%
    rownames_to_column("FID_IID") %>%
    separate(FID_IID, c("FID", "IID"), sep = sep)
}


#' Take a matrix product of data matrix (X_df) and coefficient matrix (beta_df)
#'
#' Given a data matrix (X_df) and a coefficient matrix (beta_df),
#' compute the matrix product.
#' You can specify the set of variables (variables) to focus
#' if you don't want to use all the data in the beta_df.
#'
#' @param X_df A data frame containing the dataset (X)
#' @param beta_df A string representing a regression formula
#' @param variables A subset of variables you'd like to take the product on
#' @param beta_estimate_cols The set of column names of the beta_df that contains the values of the coefficients.
#' @param beta_variable_col The column name of the beta_df that contains the variable name.
#' @return a data frame containing the results of XB
#' @examples
#' \dontrun{
#' compute_matrix_product(data_df, glm_fit_df, c('age', 'sex'), c('estimate'), 'variable')
#' }
#'
#' @export
compute_matrix_product <- function(
    X_df, beta_df, variables=NULL,
    beta_estimate_cols=c('estimate'),
    beta_variable_col='variable'
){
    if(is.null(variables)){
        # if set of variables are not specified,
        # we use all of the provided variables in beta_df
        beta_df %>%
        pull(all_of(beta_variable_col)) -> variables
    }

    # take the matrix product XB
    as.matrix(
        # prepater matrix X
        X_df %>% select(all_of(variables))
    ) %*% as.matrix(
        # prepare matrix B
        beta_df %>%
        rename(!!'variable_' := all_of(beta_variable_col)) %>%
        filter(variable_ %in% variables) %>%
        column_to_rownames('variable_') %>%
        select(all_of(beta_estimate_cols))
    ) %>%
    as.data.frame()
}

