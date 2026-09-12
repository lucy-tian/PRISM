library(r2redux)
library(devtools)
library(dplyr)
args<-commandArgs(TRUE)

phenotypes = readLines("/home/lucytian/data/1_Single_Cell_PRS/2_cV2F/pheno_tissue.txt")
pop = args[1]
for (p in phenotypes) {
    filename <- paste0('by_trait_', p, '_', pop, ".tsv")
    data <- read.delim(filename, sep="\t", header=TRUE)
    output=r2_diff(data,c(2),c(1),nrow(data))
    df <- data.frame(output)
    df$Phenotype <- c(p)
    write.csv(df, paste0('by_trait_', p, '_', pop, "_output.csv"), row.names = FALSE)
}