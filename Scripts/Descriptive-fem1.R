rm(list=ls())
setwd("~/Dropbox/Self-sperm-code-checking/")
library(Biobase)
library(RColorBrewer)

load("Data/eset_vst_fem1.RData")
eset <- eset.fem1

pca <- prcomp(t(exprs(eset)), scale = FALSE)
s <- summary(pca)$importance[, 1:4]

## Identify genes driving each PC:
PCgenes <- as.data.frame(pca$rotation)
write.table(PCgenes, sep = "\t", quote = FALSE, file = "Results/PCgenes_fem1.txt")

##############################################################
## color coded by age/genotype and shapes for mated/unmated ##
##############################################################
legend <- unique(pData(eset)[, c("Genotype", "Age", "col", "bg", "shp", "Mated")])

pdf("Results/pca_vst_fem1.pdf",
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
legend(x=50, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste(legend$Genotype, "day", legend$Age, legend$Mated),
       bty = "n")
dev.off()

pdf("Results/pca_vst_PC3-4_fem1.pdf",
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
legend(x=30, y=0,
       col = legend$col,
       pt.bg = legend$bg,
       pch = legend$shp,
       pt.cex = 1.3,
       legend = paste(legend$Genotype, "day", legend$Age, legend$Mated),
       bty = "n")
dev.off()
