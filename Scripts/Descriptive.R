###########################################################
## Create a PCA and other descriptive graphs of the data ##
###########################################################

rm(list=ls())
setwd("~/Dropbox/Code-For-Reviewers/")
library(Biobase)
library(RColorBrewer)

## Load VST normalized data from Normalize.R:
load("Data/eset_vst.RData")

## Perform principal component analysis:
pca <- prcomp(t(exprs(eset)), scale = TRUE)
s <- summary(pca)$importance[, 1:4]

## Identify genes driving each PC:
PCgenes <- as.data.frame(pca$rotation)
write.table(PCgenes, sep = "\t", quote = FALSE, file = "Data/PCgenes_all.txt")

#####################################################################################################
## Create a graph of the PCA data that is color coded by age/genotype and shapes for mated/unmated ##
#####################################################################################################
legend <- unique(pData(eset)[, c("Genotype", "Age", "col", "bg", "Mated", "shp")])

pdf("Results/pca_vst.pdf",
    width = 6.0, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 15.1),
    xpd = TRUE)
plot(pca$x[, "PC1"],
     pca$x[, "PC2"],
     pch = pData(eset)$shp,
     cex = 1.3,
     col = pData(eset)$col,
     bg = pData(eset)$bg,
     las = 1,
     xlab=paste("PC1 (", round(100*s[2,1], digits = 1), "%)", sep = ""),
     ylab=paste("PC2 (", round(100*s[2,2], digits = 1),  "%)", sep = ""))
legend(x=150, y=50,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste(legend$Genotype, "day", legend$Age, legend$Mated),
       bty = "n")
dev.off()

pdf("Results/pca_vst_PC2-3.pdf",
    width = 6.0, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 15.1),
    xpd = TRUE)
plot(pca$x[, "PC2"],
     pca$x[, "PC3"],
     pch = pData(eset)$shp,
     cex = 1.3,
     col = pData(eset)$col,
     bg = pData(eset)$bg,
     las = 1,
     xlab=paste("PC2 (", round(100*s[2,2], digits = 1), "%)", sep = ""),
     ylab=paste("PC3 (", round(100*s[2,3], digits = 1),  "%)", sep = ""))
legend(x=100, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste(legend$Genotype, "day", legend$Age, legend$Mated),
       bty = "n")
dev.off()

pdf("Results/pca_vst_PC3-4.pdf",
    width = 6.0, height = 3.2)
par(mar = c(4.1, 4.1, 1.1, 15.1),
    xpd = TRUE)
plot(pca$x[, "PC3"],
     pca$x[, "PC4"],
     pch = pData(eset)$shp,
     cex = 1.3,
     col = pData(eset)$col,
     bg = pData(eset)$bg,
     las = 1,
     xlab=paste("PC3 (", round(100*s[2,3], digits = 1), "%)", sep = ""),
     ylab=paste("PC4 (", round(100*s[2,4], digits = 1),  "%)", sep = ""))
legend(x=100, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste(legend$Genotype, "day", legend$Age, legend$Mated),
       bty = "n")
dev.off()

#######################################################################################################
## Hierarchical clustering, based on euclidean distance of the samples, clustering method: complete: ##
#######################################################################################################
pdf("Results/dendrogram_euclidean_complete_vst.pdf",
    width = 10, height = 5)
par(mar=c(12.1, 4.1, 1.1, 1.1))
dd <- as.dendrogram(
  hclust(
    dist(t(exprs(eset)), method = "euclidean")
  )
)
plot(dd)
dev.off()
