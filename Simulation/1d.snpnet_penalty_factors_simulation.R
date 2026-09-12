fullargs <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

script_name <- normalizePath(
    sub("--file=", "", fullargs[grep("--file=", fullargs)])
)
script_dir <- dirname(script_name)
#script_dir <- getwd()

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))


####################################################################
#source(file.path(dirname(script_dir), "paths.sh"))
source(file.path("paths.sh"))
####################################################################


####################################################################
# output files
####################################################################

# Store args[1] and args[2] in variables
arg1 <- as.numeric(args[1])

#pfactor_f <- file.path(
#    p_factors,
#    paste("ukb_genoHM3_", arg1, "_cCRE.pfactor.rds", sep="")
#)

pfactor_f <- file.path(
    "penalty_files",
    paste0("set_P_q", as.integer(arg1 * 100), ".rds")
)


#pfactor_f <- file.path(
#    p_factors,
#    paste(args[1], "_ukb_geno.snpnet.pfactor.rds", sep="")
#)

####################################################################
# input files
####################################################################

annot_f <- paste0("set_P_q", as.integer(arg1 * 100), ".tsv")


####################################################################
# main
####################################################################

annot_f %>%
fread() %>%
rename_with(~str_replace(., "#", ""), starts_with("#")) -> annot_df


annot_df %>%
colnames


annot_w_df <- annot_df %>%
    mutate(
        snpnet_w = if_else(score == 1, 0.7, 1)
    )

#mutate(
#    snpnet_w = case_when(
#      (`CA-H3K4me3` == 1 | `CA-CTCF` == 1 | `CA-TF` == 1 ) ~ 0.7,
#      (`PLS` == 1 | `pELS` == 1 | `dELS` == 1 | `CA-only` == 1) ~ 0.8,
#      (`Low-DNase` == 1) ~ 1.2,
#      TRUE ~ 1
#    )
#  ) -> annot_w_df

annot_w_df %>%
count(snpnet_w)


annot_w_df %>%
mutate(ID_ALT = paste(ID, ALT, sep='_')) %>%
select(ID_ALT, snpnet_w) %>%
deframe() %>%
saveRDS(file = pfactor_f)
