fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

save_d <- args[1]

library(dplyr)
library(data.table)
library(tidyr)

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

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

load(file.path(save_d, "snpnet.RData"))

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