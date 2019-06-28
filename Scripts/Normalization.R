rm(list=ls())
setwd("~/Dropbox/Self-sperm-code-checking/")
library(Biobase)
library(DESeq2)

##########################################
##### Normalization of all data sets #####
##########################################
## Load the eset created using Import.R:
load("Data/eset_counts.RData")

## Create a column that describes the way to group the samples for differential gene expression analysis
## (For example, group the replicates together):
pData(eset)$group <- factor(paste(pData(eset)$Genotype, pData(eset)$Mated, pData(eset)$Age, sep = "_"))

## Use DEseq2 to measure differential gene expression:
dds <- DESeqDataSetFromMatrix(countData = exprs(eset),
                              design = ~ group,
                              colData = pData(eset))
dds.deseq <- DESeq(dds)
save(dds.deseq,
     file = "Data/dds_deseq.RData")


## Estimate data dispersion:
pdf("Data/dispersion_estimate.pdf")
plotDispEsts(dds.deseq)
dev.off()

## Normalize data using Variance Stabilizing Transformation and save (this is what will be used for the PCA):
vst <- varianceStabilizingTransformation(dds.deseq,
                                         blind = TRUE)
exprs(eset) <- assay(vst)
save(eset,
     file = "Data/eset_vst.RData")
write.table(exprs(eset),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst.txt")


#####################################################
##### Normalization of N2 data sets (all genes) #####
#####################################################
rm(list=ls())
setwd("~/Documents/RNA-seq/2018-05-MID-N2fem-1/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

## Remove the fem-1 samples from the eset
eset.N2 <- eset[, c(17:32)]

## Normalize data from only the N2 samples:
pData(eset.N2)$group <- factor(paste(pData(eset.N2)$Genotype, pData(eset.N2)$Age, sep = "_"))
dds.N2 <- DESeqDataSetFromMatrix(countData = exprs(eset.N2),
                                      design = ~ group,
                                      colData = pData(eset.N2))
dds.deseq.N2 <- DESeq(dds.N2)
save(dds.deseq.N2,
     file = "Data/dds_deseq_N2.RData")


pdf("Data/dispersion_estimate_N2.pdf")
plotDispEsts(dds.deseq.N2)
dev.off()

vst.N2 <- varianceStabilizingTransformation(dds.deseq.N2,
                                                 blind = TRUE)
exprs(eset.N2) <- assay(vst.N2)
save(eset.N2,
     file = "Data/eset_vst_N2.RData")
write.table(exprs(eset.N2),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst_N2.txt")

########################################################
##### Normalization of fem-1 data sets (all genes) #####
########################################################
rm(list=ls())
setwd("~/Documents/RNA-seq/2018-05-MID-N2fem-1/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

## Remove the N2 samples from the eset
eset.fem1 <- eset[, c(1:16)]

## Normalize data from only the fem-1 samples:
pData(eset.fem1)$group <- factor(paste(pData(eset.fem1)$Genotype, pData(eset.fem1)$Age, sep = "_"))
dds.fem1 <- DESeqDataSetFromMatrix(countData = exprs(eset.fem1),
                                      design = ~ group,
                                      colData = pData(eset.fem1))
dds.deseq.fem1 <- DESeq(dds.fem1)
save(dds.deseq.fem1,
     file = "Data/dds_deseq_fem1.RData")


pdf("Data/dispersion_estimate_fem1.pdf")
plotDispEsts(dds.deseq.fem1)
dev.off()

vst.fem1 <- varianceStabilizingTransformation(dds.deseq.fem1,
                                                 blind = TRUE)
exprs(eset.fem1) <- assay(vst.fem1)
save(eset.fem1,
     file = "Data/eset_vst_fem1.RData")
write.table(exprs(eset.fem1),
            col.names = NA,
            sep = "\t",
            quote = FALSE,
            file = "Data/vst_fem1.txt")

