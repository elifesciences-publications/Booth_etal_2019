#!/usr/bin/Rscript

#This script will read in the differentially expressed gene lists and output genes that pass an FDR threshold of 0.05 that are either
#upregulated or downregulated per condition



rm(list=ls())

DEG_lists <- list.files(path="~/Dropbox/lauren_motif_enrichment_6_14_19_testing/RNA-seq_Data/DEG_files")

setwd("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/RNA-seq_Data")

for (DEGs in DEG_lists) {
  df <- read.table(DEGs)
  df <- df[df$padj <= 0.05,]
  df_pos <- df[df$log2FoldChange > 0,]
  df_neg <- df[df$log2FoldChange < 0,]
  
  DEGs <- gsub(".txt",'',DEGs)
  conditions <- strsplit(DEGs,split="v")
  write.table(df_pos, file = paste("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/RNA-seq_Data/processed_DEG_files/",conditions[[1]][1],"_up_",DEGs,".txt", sep=""), quote = FALSE, sep="\t")
  write.table(df_neg, file = paste("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/RNA-seq_Data/processed_DEG_files/",conditions[[1]][2],"_up_",DEGs,".txt", sep=""), quote = FALSE, sep="\t")
}
