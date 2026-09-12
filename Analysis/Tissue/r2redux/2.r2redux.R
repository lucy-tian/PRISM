library(r2redux)
library(devtools)
library(dplyr)
args<-commandArgs(TRUE)

phenotypes = readLines("/home/lucytian/data/1_Single_Cell_PRS/2_cV2F/pheno_tissue.txt")
pop = args[1]
#phenotypes = readLines("pheno.txt")
tissues <- c('BLOOD', 'LIVER', 'LUNG', 'KIDNEY', 'all')
for (p in phenotypes) {
    filename <- paste0(p, '_', pop, ".tsv")
    data <- read.delim(filename, sep="\t", header=TRUE)
    df_list <- list()
    for (i in 1:5) {
        output=r2_diff(data,c(i+1),c(1),nrow(data))
        df <- data.frame(output)
        df$Phenotype <- c(p)
        df$Tissue <- c(tissues[i])
        df_list <- append(df_list, list(df))
    }
    df_all <- bind_rows(df_list)
    write.csv(df_all, paste0(p, '_', pop, "_output.csv"), row.names = FALSE)
}