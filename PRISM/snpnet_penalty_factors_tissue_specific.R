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
arg1 <- args[1]
arg2 <- as.numeric(args[2])
arg3 <- args[3]

#pfactor_f <- file.path(
#    p_factors,
#    paste("ukb_genoHM3_", arg1, "_cCRE.pfactor.rds", sep="")
#)

pfactor_f <- file.path(
    p_factors,
    paste(arg1, "_cV2F.pfactor.rds", sep="")
)

dir.create(dirname(pfactor_f), recursive = TRUE, showWarnings = FALSE)

#pfactor_f <- file.path(
#    p_factors,
#    paste(args[1], "_ukb_geno.snpnet.pfactor.rds", sep="")
#)

####################################################################
# input files
####################################################################

annot_f <- arg3


####################################################################
# main
####################################################################

annot_f %>%
fread() %>%
rename_with(~str_replace(., "#", ""), starts_with("#")) -> annot_df


annot_df %>%
colnames


annot_df %>%
    mutate(
    snpnet_w = case_when(
        (
            (
                (Csq_group == "PTVs") |
                (ClinVar == "pathogenic")
            )
        ) ~ .5,
        (
            (
                (geno_source == "hla") |
                (Csq_group == "PAVs") |
                (ClinVar == "likely_pathogenic")
            )
        ) ~ .6,
                ( 
            !!sym(paste('cV2F', arg1, sep="_")) >= arg2
        ) ~ 0.7,
        (
            (!in_LDSC_hm3)
        ) ~ 1.2,
        TRUE ~ 1
    )
) -> annot_w_df

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
#write.csv(file = str_replace(pfactor_f, ".rds$", ".csv"), row.names = TRUE)
#saveRDS(file = pfactor_f)

new_col_name <- paste('cV2F', arg1, 'prioritize', sep = "_")
existing_col_name <- paste('cV2F', arg1, sep = "_")

annot_w_df_c <- annot_w_df %>%
  mutate(!!new_col_name := !!sym(existing_col_name) >= arg2)


annot_w_df_c %>%
count(
    snpnet_w,
    geno_source,
    Csq_group,
    ClinVar,
    in_LDSC_hm3,
    !!sym(new_col_name)
) -> annot_w_count_df


annot_w_count_df %>% 
rename("#snpnet_w" = "snpnet_w") %>%
fwrite(
    str_replace(pfactor_f, ".rds$", ".counts.tsv.gz"),
    sep = "\t", na = "NA", quote = F
)

