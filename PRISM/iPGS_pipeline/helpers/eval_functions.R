#' compute the difference between the full model and covariate-only model
#'
compute_delta_eval <- function(eval_df) {
    eval_df %>%
    filter(model == "full") %>%
    select(-model, -predictors, -l_eval, -u_eval, -P) %>%
    left_join(
        eval_df %>%
        filter(model == "covars") %>%
        select(
            refit_str, method, response, population, split, metric, `eval`
        ) %>%
        rename("eval_covars" = "eval"),
        by = c(
            "refit_str", "method", "response",
            "population", "split", "metric"
        )
    ) %>%
    mutate(
        `eval` = `eval` - eval_covars,
        model = "delta"
    ) %>%
    select(-eval_covars) -> eval_delta_df

    bind_rows(
        eval_df,
        eval_delta_df
    )
}



#' evaluate R2 or AUROC across different population/split
#'
#' @param df a data frame containing individual-level data
#' @param score_col the column name of the predictor variable
#' @param phenotype the column name of the response variable
#' @param family gaussian or binomial
#' @param pop population
#'
#' This function takes in a data frame containing individual-level data,
#' the column name of the predictor variable, the column name of the
#' response variable, the type of model (gaussian or binomial),
#' the population of interest.
score_eval_wrapper <- function(
    df, score_col, phenotype, family, pop
) {

    # Check that the family parameter is valid
    stopifnot(family %in% c("gaussian", "binomial"))

    pops <-
        df %>%
        filter(split == "test") %>%
        pull(population) %>%
        unique()
    # Check that the pop parameter is valid
    stopifnot((pop == "all") || (pop %in% pops))

    if (pop == "all") {
        bind_rows(
            # Evaluate on train and val splits together
            tryCatch({
                df %>%
                filter(split %in% c("train", "val")) %>%
                eval_CI(phenotype, score_col, family) %>%
                mutate(split = "train_val")
            }, error = function(e) {message(e)}),

            # Evaluate on train, val, and test splits separately
            c("train", "val", "test") %>%
            lapply(function(s) {
                tryCatch({
                    df %>%
                    filter(split == s) %>%
                    eval_CI(phenotype, score_col, family) %>%
                    mutate(split = s)
                }, error = function(e) {message(e)})
            }) %>%
            bind_rows(),

            # Evaluate on different populations in the test split
            lapply(pops, function(p) {
                tryCatch({
                    df %>%
                    filter(split == "test", population == p) %>%
                    eval_CI(phenotype, score_col, family) %>%
                    mutate(split = "test", population = p)
                }, error = function(e) {message(e)})
            }) %>%
            bind_rows(),
            
            lapply(pops, function(p) {
                tryCatch({
                    df %>%
                    filter(split == "val", population == p) %>%
                    eval_CI(phenotype, score_col, family) %>%
                    mutate(split = "val", population = p)
                }, error = function(e) {message(e)})
            }) %>%
            bind_rows(),            
            
        ) -> eval_df

    } else {
        bind_rows(
            tryCatch({
                df %>%
                filter(split %in% c("train", "val"), population == pop) %>%
                eval_CI(phenotype, score_col, family) %>%
                mutate(split = "train_val")
            }, error = function(e) {message(e)}),

            c("train", "val", "test") %>%
            lapply(function(s) {
                tryCatch({
                    df %>%
                    filter(split == s, population == pop) %>%
                    eval_CI(phenotype, score_col, family) %>%
                    mutate(split = s)
                }, error = function(e) {message(e)})
            }) %>%
            bind_rows()
        ) %>%
        mutate(
            population = pop
        ) -> eval_df
    }

    return(eval_df)
}

get_plot_label <- function(GBE_ID, GBE_names = setNames(NULL, NULL)) {
    if (GBE_ID %in% names(GBE_names))
        sprintf("%s: %s", GBE_ID, GBE_names[[GBE_ID]])
    else
        GBE_ID
}

####################################################################
# A copy of psychometric::CI.Rsq() function written by Thomas D. Fletcher
#
# https://rdrr.io/cran/psychometric/man/CI.Rsq.html
# https://rdrr.io/cran/psychometric/src/R/CI.Rsq.R
#
# As this is the only function we use from psychometric package (GPL2),
# we place a copy here for now.
#
# 2021.11.29
####################################################################
"CI.Rsq" <-
function(rsq, n, k, level=.95)
 {
noma <- 1-level
sersq <- sqrt((4*rsq*(1-rsq)^2*(n-k-1)^2)/((n^2-1)*(n+3)))
zs <- - qnorm(noma/2)
mez <- zs*sersq
lcl <- rsq - mez
ucl <- rsq + mez
mat <- data.frame(Rsq = rsq, SErsq = sersq, LCL = lcl, UCL = ucl)
return(mat)
}


####################################################################
# Plot functions
####################################################################


plot_eval_CI <- function(
    eval_df, family, phenotype, GBE_names = setNames(NULL, NULL)
) {
    population_levels_all <- c(
        "all",
        eval_df %>%
            filter(
                split == "test",
                !is.na(population)
            ) %>%
            pull(population) %>%
            unique()
    )

    eval_df %>%
    filter(split == "test") %>%
    replace_na(list(population = "all")) -> eval_filtered_df

    eval_filtered_df %>%
    select(all_of(c("population"))) %>%
    unique() %>%
    pull() -> evaluated_populations

    eval_filtered_df %>%
    mutate(
        population = factor(
            population,
            levels = rev(population_levels_all[
                population_levels_all %in% c(evaluated_populations)
            ])
        )
    ) %>%
    ggplot(aes(x = `eval`, y = population, color = model)) -> plot_eval

    if (family == "binomial") {
        plot_eval <- plot_eval +
        geom_vline(xintercept = .5, color = "gray") +
        labs(x = "AUROC")
    }

    if (family == "gaussian") {
        plot_eval <- plot_eval +
        geom_vline(xintercept = 0, color = "gray") +
        labs(x = latex2exp::TeX("$\\textit{R}^2$"))
    }

    plot_eval +
    geom_errorbarh(
        aes(xmin = l_eval, xmax = u_eval),
        position = position_dodge(width = 0.6)
    ) +
    geom_point(position = position_dodge(width = 0.6)) +
    theme_bw(base_size = 16) +
    labs(
        title = get_plot_label(phenotype, GBE_names),
        y = "Population groups in UK Biobank"
    ) -> plot_eval

    plot_eval
}


#' generate a 2d heatmap plot comparing the PRS vs phenotype.
#' Z-score of PRS is ued when geno_z == TRUE
#'
#' @param plot_df PRS and phenotype columns for (typically test) set of individuals
#' @param plot_bin2d_x size of the x-axis grid in 2d heatmap
#' @param plot_bin2d_y size of the y-axis grid in 2d heatmap
#' @param geno_score_col column name of the score in plot_df
#' @param phe_col column name of the phenotype data in plot_df
#' @param geno_z whether to take the Z-score transformation of the score
#'
plot_PRS_vs_phe <- function(
    plot_df, plot_bin2d_x=NULL, plot_bin2d_y=NULL,
    geno_score_col="geno_score", phe_col="phe", geno_z=TRUE
) {
    # rename the columns
    plot_df %>%
    rename(!!"phe" := all_of(phe_col)) %>%
    rename(!!"geno_score" := all_of(geno_score_col)) -> plot_dff

    if (geno_z) {
        plot_dff %>% mutate(
            geno_score_z = (geno_score - mean(geno_score)) / sd(geno_score)
        ) -> plot_dff
        if (is.null(plot_bin2d_x)) plot_bin2d_x <- 0.05
    }else{
        plot_dff %>%
        mutate(geno_score_z = geno_score) -> plot_dff
        if (is.null(plot_bin2d_x)) {
            # compute the bin size for x-axis
            plot_bin2d_x <- diff(quantile(plot_dff$geno_score, c(.4, .6))) / 4
        }
    }

    if (is.null(plot_bin2d_y)) {
        # compute the bin size for y-axis
        plot_bin2d_y <- diff(quantile(plot_dff$phe, c(.4, .6))) / 4
    }

    # plot the bin2d plot for the 99.9% coverage
    plot_dff %>%
    filter(
        # 99.9% coverage
        quantile(plot_dff$phe, .0005) < phe,
        quantile(plot_dff$phe, .9995) > phe
    ) %>%
    ggplot(aes(x = geno_score_z, y = phe)) +
    geom_bin2d(binwidth = c(plot_bin2d_x, plot_bin2d_y)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw(base_size = 16) +
    labs(x = ifelse(geno_z, "snpnet PGS (Z-score)", "snpnet PGS"))
}


#' generate Violin plot comparing the distribution of PRSs
#'  (Z-score if geno_z == TRUE) stratified by case/control status.
#' plot_df should contain PRS and phenotype columns for
#' (typically test) set of individuals.
#'  The column names can be specified in geno_score_col and phe_col
plot_PRS_binomial <- function(
    plot_df, geno_score_col="geno_score",
    phe_col="phe", geno_z=TRUE
) {
    # rename the columns
    plot_df %>%
    rename(!!"phe" := all_of(phe_col)) %>%
    rename(!!"geno_score" := all_of(geno_score_col)) -> plot_dff

    if (geno_z) {
        plot_dff %>%
        mutate(
            geno_score_z = (geno_score - mean(geno_score)) / sd(geno_score)
        ) -> plot_dff
    }else{
        plot_dff %>%
        mutate(geno_score_z = geno_score) -> plot_dff
    }

    plot_dff %>%
    left_join(
        data.frame(
            phe_str = c("Control", "Case"),
            phe = c(0, 1),
            stringsAsFactors = F
        ),
        by = "phe"
    ) %>%
    ggplot(aes(
        x = reorder(as.factor(phe_str), phe),
        y = geno_score_z,
        color = as.factor(phe)
    )) +
    geom_violin() +
    geom_boxplot(outlier.size = 0, outlier.stroke = 0, width = 0.2) +
    stat_summary(
        fun = mean, geom = "errorbar",
        aes(ymax = ..y.., ymin = ..y..),
        width = 1.1, linetype = "dashed"
    ) +
    theme_bw(base_size = 16) +
    theme(legend.position = "none") +
    labs(
        x = "phenotype",
        y = ifelse(geno_z, "snpnet PGS (Z-score)", "snpnet PGS")
    )
}


plot_PRS_bin_vs_phe <- function(summary_plot_df, horizontal_line) {
    summary_plot_df %>%
    mutate(
        x_ticks_labels = if_else(
            endsWith(bin_str, "100%"),
            paste0("[", bin_str, "]"),
            paste0("[", bin_str, ")")
        )
    ) %>%
    ggplot(aes(x = reorder(x_ticks_labels, u_bin), y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = l_err, ymax = u_err)) +
    geom_hline(yintercept = horizontal_line, color = "gray") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    labs(x = "snpnet PGS percentile")
}


plot_PRS_bin_vs_OR <- function(summary_plot_df) {
    summary_plot_df %>%
    mutate(mean = OR, l_err = l_OR, u_err = u_OR) %>%
    plot_PRS_bin_vs_phe(1) +
    labs(y = "Odds ratio [SE]")
}

plot_PRS_bin_vs_phe_multi_methods <-
function(summary_plot_df, horizontal_line) {
    summary_plot_df %>%
    mutate(
        x_ticks_labels = if_else(
            endsWith(bin_str, "100%"),
            paste0("[", bin_str, "]"),
            paste0("[", bin_str, ")")
        )
    ) %>%
    ggplot(aes(
        x = reorder(x_ticks_labels, u_bin),
        y = mean,
        color = method,
        shape = method
    )) +
    geom_point(
        position = position_dodge(width = 0.7)
    ) +
    geom_errorbar(
        aes(ymin = l_err, ymax = u_err),
        position = position_dodge(width = 0.7)
    ) +
    geom_hline(yintercept = horizontal_line, color = "gray") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    labs(x = "snpnet PGS percentile")
}


plot_PRS_bin_vs_OR_multi_methods <-
function(summary_plot_df) {
    summary_plot_df %>%
    mutate(mean = OR, l_err = l_OR, u_err = u_OR) %>%
    plot_PRS_bin_vs_phe_multi_methods(1) +
    labs(y = "Odds ratio [SE]")
}



####################################################################
# Quantile plot functions (compute (mean and sd) or odds ratio )
####################################################################

#' @importFrom magrittr %>%
#'
compute_mean <- function(df, percentile_col, phe_col, l_bin, u_bin) {
    # Compute the mean and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    stratified_df <- df %>%
    rename(
        !!"Percentile" := all_of(percentile_col),
        !!"phe" := all_of(phe_col)
    ) %>%
    filter(l_bin < Percentile, Percentile <= u_bin)

    n     <- stratified_df %>%
        nrow()
    mean  <- stratified_df %>%
        select(phe) %>%
        pull() %>%
        mean()
    sd    <- stratified_df %>%
        select(phe) %>%
        pull() %>%
        sd()
    std_e <- sd / sqrt(n)
    l_err <- mean - std_e
    u_err <- mean + std_e

    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        mean   = mean,
        std_err = std_e,
        l_err = l_err,
        u_err = u_err,
        mean_str = sprintf("%.3f (%.3f,%.3f)", mean, l_err, u_err),
        bin_str = paste0(100 * l_bin, "% - ", 100 * u_bin, "%"),
        stringsAsFactors = F
    ) %>%
    mutate(mean_str = as.character(mean_str))
}


#' a helper function for compute_OR.
#' This provides the counts of the descrete phenotype value (phe_col)
#' for the specified bin (l_bin, u_bin], based on the percentile of PRS (percentile_col)
#'
#' @importFrom magrittr %>%
#'
filter_by_percentile_and_count_phe <- function(
    df, percentile_col, phe_col, l_bin, u_bin
) {
    df %>%
    rename(
        !!"Percentile" := all_of(percentile_col),
        !!"phe" := all_of(phe_col)
    ) %>%
    filter(l_bin < Percentile, Percentile <= u_bin) %>%
    count(phe) %>%
    # To cope with sparse bins where case or control counts are zero,
    # we add the following dummy counts (zeros)
    bind_rows(data.frame(
        phe = c(0, 1),
        n = as.integer(c(0, 0))
    )) %>%
    group_by(phe) %>%
    summarise(n = sum(n)) %>%
    ungroup
}


#' Compute the OR and sd of the trait value (phe_col), based on the
#' binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
#' The odds ratio is defined against the "middle" of the PRS distribution, and
#' we assume to have the phenotype counts in that bin (cnt_middle)
#'
#' @importFrom magrittr %>%
#'
compute_OR <- function(df, percentile_col, phe_col, l_bin, u_bin, cnt_middle) {
    cnt_tbl <- df %>%
    filter_by_percentile_and_count_phe(
        percentile_col, phe_col, l_bin, u_bin
    ) %>%
    bind_rows(cnt_middle %>% mutate(n = 0) %>%
    mutate(n = as.integer(n))) %>%
    group_by(phe) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    inner_join(cnt_middle, by = "phe") %>%
    gather(bin, cnt, -phe) %>%
    arrange(-phe, bin)

    cnt_res <- cnt_tbl %>%
        mutate(cnt = as.numeric(cnt)) %>%
        select(cnt) %>%
        pull()
    names(cnt_res) <- c("n_TP", "n_FN", "n_FP", "n_TN")

    OR <- (cnt_res[["n_TP"]] * cnt_res[["n_TN"]]) / (cnt_res[["n_FP"]] * cnt_res[["n_FN"]])
    LOR <- log(OR)
    se_LOR <- cnt_tbl %>%
        select(cnt) %>%
        pull() %>%
        lapply(function(x){1/x}) %>%
        reduce(function(x, y){x+y}) %>%
        sqrt()
    l_OR <- exp(LOR - 1.96 * se_LOR)
    u_OR <- exp(LOR + 1.96 * se_LOR)

    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        n_TP = cnt_res[["n_TP"]],
        n_FN = cnt_res[["n_FN"]],
        n_FP = cnt_res[["n_FP"]],
        n_TN = cnt_res[["n_TN"]],
        OR   = OR,
        SE_LOR = se_LOR,
        l_OR = l_OR,
        u_OR = u_OR,
        OR_str = sprintf("%.3f (%.3f,%.3f)", OR, l_OR, u_OR),
        bin_str = paste0(100 * l_bin, "% - ", 100 * u_bin, "%"),
        stringsAsFactors = F
    )
}


#' @importFrom magrittr %>%
#'
compute_summary_OR_df <- function(df, percentile_col, phe_col, bins=NULL) {
    if (is.null(bins)) bins <- ((0:10)/10)
    cnt_middle <- df %>%
    filter_by_percentile_and_count_phe(percentile_col, phe_col, 0.4, 0.6) %>%
    rename("n_40_60" = "n")

    1:(length(bins)-1) %>%
    lapply(function(i) {
        compute_OR(df, percentile_col, phe_col, bins[i], bins[i+1], cnt_middle)
    }) %>%
    bind_rows()
}



#' @importFrom magrittr %>%
#'
compute_summary_mean_df <- function(df, percentile_col, phe_col, bins=NULL) {
    if (is.null(bins)) bins <- c(0, .0005, .01, (1:19)/20, .99, .9995, 1)
    1:(length(bins)-1) %>%
    lapply(function(i) {
        compute_mean(df, percentile_col, phe_col, bins[i], bins[i+1])
    }) %>%
    bind_rows()
}


compute_summary_df <- function(
    df, percentile_col, phe_col, bins=NULL, family="gaussian"
) {
    if (family == "gaussian") {
        compute_summary_mean_df(df, percentile_col, phe_col, bins)
    }else if (family == "binomial") {
        compute_summary_OR_df(df, percentile_col, phe_col, bins)
    }else{
        stop(sprintf("%s family is not supported!", family))
    }
}



####################################################################
# Predictive performance evaluation functions (compute R2 or AUC)
####################################################################

#' Get r-squared value from glm fit object
#'
#' @param glm_fit A fit object from glm()
#' @return  The r-squared value
#'
glm_fit_to_R2 <- function(glm_fit) {
    with(summary(glm_fit), 1 - deviance / null.deviance)
}


#' Compose a regression formula give a response variable and the predictors
#'
#'  response ~ 1 + predictor[1] + predictor[2] + ... + predictor[n]
#'
#' The function quotes variable names using quote_char
#'
#' @param response The response variable
#' @param predictors The predictor variables
#' @param quote_char A character for quotation of terms
#' @return  A string representing the regressino formula
#' @examples
#' compose_regression_formula_str('HC326', c('covar', 'PRS.score-1'))
#' compose_regression_formula_str('HC326', c('covar', 'PRS.score-1'), quote_char='')
#'
#' @export
compose_regression_formula_str <- function(
    response, predictors, quote_char="`"
) {
    return(sprintf(
        "%s ~ 1 + %s",
        paste0(quote_char, response, quote_char),
        paste(sapply(
            predictors,
            function(term){paste0(quote_char, term, quote_char)}
        ), collapse = " + ")
    ))
}


#' Compute AUC with confidence interval.
#'
#' Given a data frame containing a binary response variable and predictor(s),
#' compute the AUC with confidence interval and the p-value.
#' If multiple predictors are specified, we return the minimum p-value
#'
#' @param data A data frame containing the response and predictor(s)
#' @param response The name of the response variable
#' @param predictors The predictor variable(s)
#' @param level The confidence level of the interval
#' @return a data frame containing AUC (eval),
#' confidence interval (l_eval, u_eval), p-value along with
#' the information of the specified response and predictors
#' @examples
#' \dontrun{
#' eval_AUROC_CI(data_df, 'HC326', c('covar_score'))
#' }
#'
#' @importFrom magrittr %>%
#'
#' @export
eval_AUROC_CI <- function(data, response, predictors, level=.95) {
    # regression formula
    formula_str <- compose_regression_formula_str(response, predictors)

    # get the p-value
    if (all(unique(data[[response]]) %in% c(1, 2))) {
        data[[response]] <- data[[response]] - 1
    }
    
    data %>% fit_glm(formula_str, "binomial") %>% fit_to_df() %>%
    mutate(variable = str_replace_all(variable, "`", "")) %>%
    filter(variable %in% predictors) %>% pull(P) %>%
    # we extract the smallest p-values across multiple predictors for now
    min() -> P_val

    # pROC::roc does not handle the quoted column names
    # this causes issue when predictor and response variable names contain characters like '-'
    # because it is interpreted as a minus operater not a part of the variable name.
    # To handle this, we rename column names when computing AUC-ROC with pROC::roc()
    # Note the glm() function handles such variable names as long as it's quoted.
    col_rename_func <- function(col_name){str_replace_all(col_name, "[-]", ".")}

    # call pROC::roc and pROC::ci.auc to get the AUC with CI
    pROC::ci.auc(
        pROC::roc(
            formula = stats::as.formula(
                compose_regression_formula_str(response, col_rename_func(predictors), quote_char = "")
            ),
            data = (data %>% rename_with(col_rename_func)),
            direction="<", levels=c("control"=0, "case"=1)
        ), method="delong",  conf.level=level
    ) -> AUC_ci_list

    # format the results as a data frame
    data.frame(
        response = response,
        predictors = paste(predictors, collapse = "+"),
        metric = "AUROC",
        `eval` = AUC_ci_list[2],
        l_eval = AUC_ci_list[1],
        u_eval = AUC_ci_list[3],
        P = P_val,
        n_case    = data %>%
            rename(!!"phe" := all_of(response)) %>%
            filter(phe == 1) %>% nrow(),
        n_control = data %>%
            rename(!!"phe" := all_of(response)) %>%
            filter(phe == 0) %>% nrow(),
        n = nrow(data)
    )
}


#' Compute R-squared with confidence interval.
#'
#' Given a data frame containing a continupus response variable and predictor(s),
#' compute the R-squared with confidence interval and the p-value.
#' If multiple predictors are specified, we return the minimum p-value
#'
#' @param data A data frame containing the response and predictor(s)
#' @param response The name of the response variable
#' @param predictors The predictor variable(s)
#' @param level The confidence level of the interval
#' @return a data frame containing r-squared (eval),
#'   confidence interval (l_eval, u_eval), p-value along with
#'   the information of the specified response and predictors
#' @examples
#' \dontrun{
#' eval_R2_CI(data_df, 'INI50', c('covar_score'))
#' }
#'
#' @export
eval_R2_CI <- function(data, response, predictors, level=.95) {
    data %>% fit_glm(
        compose_regression_formula_str(response, predictors),
        "gaussian"
    ) -> glm_fit
    # get the p-value
    glm_fit %>% fit_to_df() %>%
    mutate(variable = str_replace_all(variable, "`", "")) %>%
    filter(variable %in% predictors) %>%
    pull(P) %>%
    # we extract the smallest p-values across multiple predictors for now
    min() -> P_val
    # compute the r-squared value
    glm_fit %>% glm_fit_to_R2() -> rsq
    # call psychometric::CI.Rsq() to compute confidence interval
    # https://rdrr.io/cran/psychometric/man/CI.Rsq.html
    # https://rdrr.io/cran/psychometric/src/R/CI.Rsq.R
    CI.Rsq(rsq, n=nrow(data), k=length(predictors), level=level) %>%
    # format the resulting table
    select(-SErsq) %>% mutate(
        metric = "R2",
        response = response,
        predictors = paste(predictors, collapse = "+"),
        P = P_val,
        n = nrow(data)
    ) %>%
    rename("eval"="Rsq", "l_eval"="LCL", "u_eval"="UCL") %>%
    select(response, predictors, metric, `eval`, l_eval, u_eval, P, n)
}


#' Compute R-squared or AUC with confidence interval.
#'
#' Given a data frame containing a response variable and predictor(s),
#' compute R-squared or AUC with confidence interval and the p-value.
#' We currently supports gaussian and binomial models.
#' If multiple predictors are specified, we return the minimum p-value
#'
#' @param data A data frame containing the response and predictor(s)
#' @param response The name of the response variable
#' @param predictors The predictor variable(s)
#' @param level The confidence level of the interval
#' @return a data frame containing r-squared or AUC (eval),
#'   confidence interval (l_eval, u_eval), p-value along with
#'   the information of the specified response and predictors
#' @examples
#' \dontrun{
#' eval_CI(data_df, 'INI50', c('covar_score'), "gaussian")
#' eval_CI(data_df, 'HC326', c('covar_score'), "binomial")
#' }
#'
#' @export
eval_CI <- function(data, response, predictors, family, level=.95) {
    stopifnot(family %in% c("gaussian", "binomial"))
    if (family == "gaussian") {
        eval_R2_CI(data, response, predictors, level)
    }else{
        eval_AUROC_CI(data, response, predictors, level)
    }
}
