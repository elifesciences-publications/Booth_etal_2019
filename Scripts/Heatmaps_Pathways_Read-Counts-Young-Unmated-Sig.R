## Creates heat maps using the EST normalized read counts

rm(list=ls())
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
setwd("~/Dropbox/Self-sperm-code-checking/")
load("Data/eset_counts.RData")
load("Data/Cele_GO_terms.RData")

######################################################
##### Normalization of young, no males data sets #####
######################################################
## Remove the mated and middle-aged samples from the eset
eset.y.unmated <- eset[, c(5:8, 21:24)]

## Normalize data from only the young unmated samples:
dds.y.unmated <- DESeqDataSetFromMatrix(countData = exprs(eset.y.unmated),
                                    design = ~ Genotype,
                                    colData = pData(eset.y.unmated))
dds.deseq.y.unmated <- DESeq(dds.y.unmated)

vst.y.unmated <- varianceStabilizingTransformation(dds.deseq.y.unmated,
                                               blind = TRUE)
exprs(eset.y.unmated) <- assay(vst.y.unmated)
save(eset.y.unmated,
     file = "Data/eset_vst_y_unmated.RData")


## Sort by genotype (N2 vs. fem-1):
eset <- eset.y.unmated[, order(pData(eset.y.unmated)$Genotype, decreasing = FALSE)]

## Rename GO lists:
gene.sets <- GO.list

##################################################
## Significant, differentially expressed genes? ##
##################################################
## Uses the DEseq2 results to make a list of significant gene expression changes
N3UvF3U <- read.delim("Results/DEseq/Groups/N3UvF3U.txt", 
                      stringsAsFactors = FALSE,
                      row.names = 1)

sig <- which(N3UvF3U[, "padj"] < 0.05)
sig.genes <- row.names(N3UvF3U)[sig]


#####################################################################
##### Heatmaps of young, unmated worm RNA-seq normalized counts #####
#####################################################################
## For each heatmap, the normalized eset data is subset for the genes in a given GO group
## Colors are normalized by row

#### GO.0006412..translation
genes <- gene.sets[["GO:0006412--translation"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      105       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0006412-translation-sig.pdf",
    width = 4, height = 2.4, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()

pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0006412-translation-sig-names.pdf",
    width = 4, height = 20, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()


## GO.0000502..proteasome.complex
genes <- gene.sets[["GO:0000502--proteasome complex"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      23       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0000502-proteasome-sig.pdf",
    width = 4, height = 1, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()

pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0000502-proteasome-sig-names.pdf",
    width = 4, height = 4, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()


#### GO.0006629..lipid.metabolic.process
genes <- gene.sets[["GO:0006629--lipid metabolic process"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##    66       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0006629-lipid-metabolic-process-sig.pdf",
    width = 4, height = 1.8, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()

pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0006629-lipid-metabolic-process-sig-names.pdf",
    width = 4, height = 10, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()


### GO.0005764..lysosome
genes <- gene.sets[["GO:0005764--lysosome"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      28       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0005764-lysosome-sig.pdf",
    width = 4, height = 1, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()

pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0005764-lysosome-sig-names.pdf",
    width = 4, height = 4, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()


### GO.0045087..innate.immune.response
genes <- gene.sets[["GO:0045087--innate immune response"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      99       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0045087-innate_immune_response-sig.pdf",
    width = 4, height = 2.5, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()

pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0045087-innate_immune_response-sig-names.pdf",
    width = 4, height = 10, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()


#### GO.0006511..ubiquitin.dependent.protein.catabolic.process
genes <- gene.sets[["GO:0006511--ubiquitin-dependent protein catabolic process"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      29       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0006511--ubiquitin-dependent-protein-catabolic-process-sig.pdf",
    width = 4, height = 1, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()

pdf("Results/Heatmaps-counts/unmated-heatmap_GO_0006511--ubiquitin-dependent-protein-catabolic-process-sig-names.pdf",
    width = 4, height = 6, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()

