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
#args[4] <- "../outputs/plink/fb-170208-snp-q30-gq30-hwe-ss0.001.traw"
#args[5] <- "fb-170208-snp-q30-gq30-hwe-ss0.001"
#setwd("/Volumes/bio-dap15.stanford.edu/zstor/2016-tiger-wgs/basic-stats/bs-test/outputs/rarefaction")
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

cutoffs <- data_frame(suffix = c("all_0.10", "all_0.15", "all_0.20", "all_0.25", "all_0.30"), cutoff = seq(0.1,0.3,0.05))

sites <- genotypes %>% select(chrom = CHR, start = POS) %>% mutate(stop = start)

filters <- foreach(i = 1:nrow(cutoffs), .combine = cbind) %dopar% {
  cutoff <- cutoffs$cutoff[i]
  filter <- apply(pop_freq[ ,2:ncol(pop_freq)], 1, function(x){all(x >= cutoff & x <= 1-cutoff)})
  write.table(sites[filter,], file = paste(args[5], "-", cutoffs$suffix[i], ".bed", sep = ""), col.names = FALSE, quote=FALSE, row.names = FALSE)
}
