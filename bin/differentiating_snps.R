#!/bin/Rscript

library(readr)
library(tidyr)
library(dplyr)
library(gplots)
library(foreach)
library(doMC)
args <- commandArgs(trailingOnly = TRUE)

## For testing
#args <- 1
#args[1] <- 4
#args[2] <- "jacksoni,sumatrae,altaica,tigris"
#args[3] <- "../sample_groups.tsv"
#args[4] <- "../outputs/plink/fb-170208-snp-q30-gq30-hwe-ss0.01.traw"
#args[5] <- "fb-170208-snp-q30-gq30-hwe-ss0.001"
###

registerDoMC(cores = args[1])
pops <- strsplit(args[2], ',')[[1]]

sample_groups <- read_tsv(args[3], col_names=c("ID", "pop"))
genotypes <- read_tsv(args[4], col_names=TRUE)
sample_groups <- sample_groups %>% filter(pop %in%  pops)
colnames(genotypes) <- gsub("_.*", "", colnames(genotypes))
genotypes_cols <- colnames(genotypes)
sample_groups <- sample_groups[sample_groups$ID %in% genotypes_cols, ]
sample_groups$idcol <- match(sample_groups$ID, genotypes_cols)

# calculate MAF per sample group
pop_freq <- foreach(pop.i = pops, .combine = cbind) %do% {
  samples.i <- sample_groups %>% filter(pop == pop.i)
  cols.i <- genotypes_cols %in% samples.i$ID
  genotypes.i <- genotypes[ , cols.i]
  genotypes.i$freq <- apply(genotypes.i,1,function(x){sum(x, na.rm = TRUE) / (sum(!is.na(x)) * 2)})
  genotypes.i <- genotypes.i %>% select(freq)
  names(genotypes.i) <- pop.i
  genotypes.i
}

# Filtering

fixed_snps <- pop_freq[apply(pop_freq, 1, function(x){all((x == 0 | x == 1) & (!is.na(x))) & ((sum(x == 1) == 1) | (sum(x==0) == 1))}), ]

bed <- genotypes %>% mutate(start = POS, stop = POS) %>% select(CHR, start, stop, SNP)

nulldrain <- foreach(pop.i = pops) %do% {
  col.i <- which(pop.i == names(fixed_snps))
  pop.i.snps <- fixed_snps[apply(fixed_snps,1,function(x){all(x[col.i] != x[-col.i])}),]
  write.table(bed[as.numeric(rownames(pop.i.snps)),], file = paste(args[5], "-diff_",pop.i, ".bed", sep = ''), col.names = FALSE, quote=FALSE, row.names = FALSE)
}

