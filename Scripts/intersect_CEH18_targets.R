#!/usr/bin/Rscript

#The purpose of this script is to output a list of ChIP-seq peaks that are associated with different lists of DEGs as well as output a unique gene list containing the DEGs that have nearby ChIP-seq peaks.

source("https://bioconductor.org/biocLite.R") #Bioconductor version 3.6 (BiocInstaller 1.28.0) R 3.4.4
biocLite("TxDb.Celegans.UCSC.ce11.refGene")
biocLite("org.Ce.eg.db")
source("https://bioconductor.org/biocLite.R")
biocLite("genomation")
source("https://bioconductor.org/biocLite.R")
biocLite("sets")

rm(list=ls())
library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Celegans.UCSC.ce11.refGene)
library(org.Ce.eg.db)
library(ChIPseeker)
library(genomation)
library(sets)
txdb <- TxDb.Celegans.UCSC.ce11.refGene

setwd("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/")


#annotate the CEH18 ChIP-seq (de-duplicated) file so the peaks are associated with nearest genes using ChIPSeeker
CEH18 <- readBed(file="~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/CEH18_ChIP_unique.bed", track.line = FALSE, remove.unusual = FALSE, zero.based = FALSE)
CEH18 <- as.data.frame(annotatePeak(CEH18, tssRegion=c(-300, 300), TxDb=txdb, annoDb="org.Ce.eg.db"))
write.table(CEH18, file = "CEH18_ChIP_annotated.bed", quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#intersecting the nearest gene CEH18 targets with the original (un-subsetted) DEG lists
setwd("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/DEGs")
DEG_list <- list.files(path="~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/DEGs/")
for (DEGs in DEG_list) {
  df <- read.table(DEGs, sep="\t", header=TRUE, row.names = NULL, quote="")
  df <- na.omit(df)
  names(df)[1] <- "GENE_NAMES"
  df_intersect <- CEH18[CEH18$SYMBOL %in% df$GENE_NAMES,]
  unique_genes <- unique(df_intersect$SYMBOL)
  
  #Manually computing Jaccard Index to determine similarity between CEH-18-associated DEGs and full DEG set
  I <- length(intersect(unique_genes,df$GENE_NAMES))
  S <- I/(length(unique_genes)+length(df$GENE_NAMES)-I)
  print(DEGs)
  print(S)
  
  DEGs <- gsub("txt",'',DEGs)
  write.table(df_intersect, file = paste("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/5_intersection_OUTPUT/",DEGs,"_intersect_CEH18_targets.bed", sep=""), quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(unique_genes, file = paste("~/Dropbox/LAUREN_motif_enrichment_6_14_19_testing/5_intersection_OUTPUT/",DEGs,"_intersect_CEH18_unique_genes.txt", sep=""), quote = FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}
