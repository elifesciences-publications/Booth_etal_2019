## Simplest model: Just define age/presence.of.males groups: GE ~ groups
## HEADS UP: If p values are reported as 0, that just means that they are actually p < 1/num.samples and should be reported like this.
## in the mean time...it is accurate to say that p<0.001

rm(list=ls())
setwd("~/Dropbox/Code-For-Reviewers/")
library(Biobase)
library(DESeq2)
load("Data/eset_counts.RData")

#For R 3.2.0
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("devtools")
#biocLite("DESeq2")
#biocLite("RcppArmadillo")


## grouping by the different variables (age, genotype, and whether they were mated)
pData(eset)$group <- factor(paste(pData(eset)$Genotype, pData(eset)$Age, pData(eset)$Mated, sep = "_"))


###################################################
## Differential expression using standard groups ##
###################################################
dds <- DESeqDataSetFromMatrix(countData = exprs(eset),
                             design = ~ group,
                             colData = pData(eset))
dds.deseq <- DESeq(dds)

## Saving the results of the DEseq2 analysis:
F7UvF3U <- results(dds.deseq, contrast = list("groupfem.1_7_unmated", "groupfem.1_3_unmated"))
F7MvF3M <- results(dds.deseq, contrast = list("groupfem.1_7_mated", "groupfem.1_3_mated"))
N7UvN3U <- results(dds.deseq, contrast = list("groupN2_7_unmated", "groupN2_3_unmated"))
N7MvN3M <- results(dds.deseq, contrast = list("groupN2_7_mated", "groupN2_3_mated"))
save(F7UvF3U, F7MvF3M, N7UvN3U, N7MvN3M,
     file = "Data/7v3_deseq.RData")

F3MvF3U <- results(dds.deseq, contrast = list("groupfem.1_3_mated", "groupfem.1_3_unmated"))
F7MvF7U <- results(dds.deseq, contrast = list("groupfem.1_7_mated", "groupfem.1_7_unmated"))
N3MvN3U <- results(dds.deseq, contrast = list("groupN2_3_mated", "groupN2_3_unmated"))
N7MvN7U <- results(dds.deseq, contrast = list("groupN2_7_mated", "groupN2_7_unmated"))
save(F3MvF3U, F7MvF7U, N3MvN3U, N7MvN7U,
     file = "Data/MvU_deseq.RData")

N3UvF3U <- results(dds.deseq, contrast = list("groupN2_3_unmated", "groupfem.1_3_unmated"))
N7UvF7U <- results(dds.deseq, contrast = list("groupN2_7_unmated", "groupfem.1_7_unmated"))
N3MvF3M <- results(dds.deseq, contrast = list("groupN2_3_mated", "groupfem.1_3_mated"))
N7MvF7M <- results(dds.deseq, contrast = list("groupN2_7_mated", "groupfem.1_7_mated"))
save(N3UvF3U, N7UvF7U, N3MvF3M, N7MvF7M,
     file = "Data/N2vfem1_deseq.RData")

all.comparisons <- list(F7UvF3U = F7UvF3U,
                        F7MvF3M = F7MvF3M,
                        N7UvN3U = N7UvN3U,
                        N7MvN3M = N7MvN3M,
                        F3MvF3U = F3MvF3U,
                        F7MvF7U = F7MvF7U, 
                        N3MvN3U = N3MvN3U,
                        N7MvN7U = N7MvN7U,
                        N3UvF3U = N3UvF3U,
                        N7UvF7U = N7UvF7U,
                        N3MvF3M = N3MvF3M,
                        N7MvF7M = N7MvF7M)

for(i in 1:length(all.comparisons)){
  comparison.i <- all.comparisons[[i]]
  file.name <- names(all.comparisons[i])
  file.path <- paste0("Results/DEseq/Groups/", file.name, ".txt")
  write.table(comparison.i,
              col.names = NA,
              sep = "\t",
              quote = FALSE,
              file = file.path)
}
