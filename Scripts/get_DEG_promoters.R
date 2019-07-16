#!/usr/bin/Rscript

#This script uses the output of "pre_processing_DEGs.R" to read in the processed DEGs (e.g. FDR-thresholded and separated by "up" or "down" per condition)
#and output the corresponding promoter regions associated with each DEG.

#NOTE: The list of promoters generated is smaller than the DEG lists since some of the gene names from the DEG lists did not appear
#in the ce11 refGene atlas.

source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Celegans.UCSC.ce11.refGene")
biocLite("org.Ce.eg.db")

rm(list=ls())
library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(org.Ce.eg.db)
library(ChIPseeker)
txdb <- TxDb.Celegans.UCSC.ce11.refGene

processed_DEG_lists <- list.files(path="~/Dropbox/lauren_motif_enrichment_6_14_19_testing/RNA-seq_Data/processed_DEG_files")

setwd("~/Dropbox/lauren_motif_enrichment_6_14_19_testing/RNA-seq_Data/processed_DEG_files")

promoters_all <- getPromoters(TxDb=txdb, upstream=300,downstream=300)
promoters_all <- as.data.frame(annotatePeak(promoters_all, tssRegion=c(-300, 300), TxDb=txdb, annoDb="org.Ce.eg.db"))
promoters_all <- na.omit(promoters_all)


for (DEGs in processed_DEG_lists) {
  df <- read.table(DEGs, sep="\t", header=TRUE, row.names = NULL)
  df <- na.omit(df)
  names(df)[1] <- "SYMBOL"
  df_prom <- promoters_all[promoters_all$SYMBOL %in% df$SYMBOL,]
  df_prom <- df_prom[,c(1:3)]
  
  DEGs <- gsub(".txt",'',DEGs)
  write.table(df_prom, file = paste("~/Dropbox/lauren_motif_enrichment_6_14_19_testing/RNA-seq_Data/promoter_DEG_lists/",DEGs,"_prom.bed", sep=""), quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

write.table(promoters_all[,c(1:3)], file = paste("~/Dropbox/lauren_motif_enrichment_6_14_19_testing/RNA-seq_Data/promoter_DEG_lists/all_promoters.bed", sep=""), quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
