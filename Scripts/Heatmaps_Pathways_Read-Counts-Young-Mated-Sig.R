## Creates heat maps using the EST normalized read counts

rm(list=ls())
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
setwd("~/Dropbox/Self-sperm-code-checking/")
load("Data/eset_counts.RData")
load("Data/Cele_GO_terms.RData")

####################################################
##### Normalization of young, +males data sets #####
####################################################
## Remove the unmated and middle-aged samples from the eset
eset.y.mated <- eset[, c(1:4, 17:20)]

## Normalize data from only the mated samples:
dds.y.mated <- DESeqDataSetFromMatrix(countData = exprs(eset.y.mated),
                                    design = ~ Genotype,
                                    colData = pData(eset.y.mated))
dds.deseq.y.mated <- DESeq(dds.y.mated)

vst.y.mated <- varianceStabilizingTransformation(dds.deseq.y.mated,
                                               blind = TRUE)
exprs(eset.y.mated) <- assay(vst.y.mated)
save(eset.y.mated,
     file = "Data/eset_vst_y_mated.RData")


## Sort by genotype (N2 vs. fem-1):
eset <- eset.y.mated[, order(pData(eset.y.mated)$Genotype, decreasing = FALSE)]

## Rename GO lists:
gene.sets <- GO.list


##################################################
## Significant, differentially expressed genes? ##
##################################################
## Uses the DEseq2 results to make a list of significant gene expression changes
N3MvF3M <- read.delim("Results/DEseq/Groups/N3MvF3M.txt", 
                      stringsAsFactors = FALSE,
                      row.names = 1)

sig <- which(N3MvF3M[, "padj"] < 0.05)
sig.genes <- row.names(N3MvF3M)[sig]



###################################################################
##### Heatmaps of young, mated worm RNA-seq normalized counts #####
###################################################################
## For each heatmap, the normalized eset data is subset for the genes in a given GO group
## The genes in the GO group are further subset to only make a heatmap with those that are DEGs (DEseq2 results, p < 0.05)
## Colors are normalized by row

#### GO.0006412..translation
genes <- gene.sets[["GO:0006412--translation"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      69       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/mated-heatmap_GO_0006412-translation-sig.pdf",
    width = 4, height = 2.4, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()


pdf("Results/Heatmaps-counts/mated-heatmap_GO_0006412-translation-sig-names.pdf",
    width = 4, height = 15, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()


### GO.0061077..chaperone.mediated.protein.folding
genes <- gene.sets[["GO:0061077--chaperone-mediated protein folding"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      5       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/mated-heatmap_GO_0061077-chaperone-mediated-protein-folding-sig.pdf",
    width = 4, height = 1, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()


pdf("Results/Heatmaps-counts/mated-heatmap_GO_0061077-chaperone-mediated-protein-folding-sig-names.pdf",
    width = 4, height = 3, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()

### GO.0040011..locomotion
genes <- gene.sets[["GO:0040011--locomotion"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      207       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/mated-heatmap_GO_0040011-locomotion-sig.pdf",
    width = 4, height = 6, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()


pdf("Results/Heatmaps-counts/mated-heatmap_GO_0040011-locomotion-sig-names.pdf",
    width = 4, height = 20, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()

### GO.0005576..extracellular.region
genes <- gene.sets[["GO:0005576--extracellular region"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      55       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/mated-heatmap_GO_0005576-extracellular-region-sig.pdf",
    width = 4, height = 4, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()


pdf("Results/Heatmaps-counts/mated-heatmap_GO_0005576-extracellular-region-sig-names.pdf",
    width = 4, height = 10, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()

### GO.0005524..ATP.binding
genes <- gene.sets[["GO:0005524--ATP binding"]]
genes <- genes[is.element(genes, rownames(eset)) & is.element(genes, sig.genes)]
eset.part <- eset[genes, ]
dim(eset.part)
## Features  Samples
##      96       8
data <- exprs(eset.part)
pdf("Results/Heatmaps-counts/mated-heatmap_GO_0005524-ATP-binding-sig.pdf",
    width = 4, height = 8, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = FALSE,
         border_color = NA)
dev.off()


pdf("Results/Heatmaps-counts/mated-heatmap_GO_0005524-ATP-binding-sig-names.pdf",
    width = 4, height = 10, onefile=FALSE)
pheatmap(data,
         scale = "row",
         color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100),
         cluster_cols = FALSE,
         show_rownames = TRUE,
         border_color = NA)
dev.off()