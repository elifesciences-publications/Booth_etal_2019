## Make expressionSet with count matrix.

setwd("~/Dropbox/Code-For-Reviewers/")
source("Scripts/make.transparent.R")
library(Biobase)

rm(list=ls())
## Create a read count matrix using each read count file from STAR:
files <- list.files("Data/Star-ReadCounts",
                    pattern = "ReadsPerGene",
                    full.names = TRUE)
samples <- gsub("Data/Star-ReadCounts/", "", files)
samples <- gsub("_ReadsPerGene.out.tab", "", samples)
system(paste("wc -l", files[1]))
## 46752 Data/Star-ReadCounts/F3M-A-GCTCATGA-AGAGGATA_S7_ReadsPerGene.out.tab
count.matrix <- matrix(nrow = 46752,
                       ncol = length(files),
                       dimnames = list(NULL, samples))

## Fill matrix with the counts for unstranded RNA-seq in the 2nd column
for(i in seq(along = files)){
  data.i <- read.table(files[i],
                       colClasses = c("NULL", "integer","NULL", "NULL"),
                       sep = "\t")[,1]
  count.matrix[, samples[i]] <- data.i
  rm(data.i)
}
genes <- read.table(files[i],
                    colClasses = c("character", "NULL","NULL", "NULL"),
                    sep = "\t")[,1]
rownames(count.matrix) <- genes


## Get Sample info for each file (make this tab-delimited file for each experiment):
pData <- read.table("Data/SampleInfo.txt",
                    sep = "\t",
                    row.names = 1,
                    header = TRUE)
pData <- pData[match(colnames(count.matrix), rownames(pData)), ]


## Assign colors and shapes to the samples (to be used when making the PCA):
source("Scripts/make.transparent.R")
pData$col[pData$Genotype == "fem-1"
          & pData$Age == 3] <- colors()[368]
pData$col[pData$Genotype == "fem-1"
          & pData$Age == 7] <- colors()[641]
pData$col[pData$Genotype == "N2"
          & pData$Age == 3] <- colors()[612]
pData$col[pData$Genotype == "N2"
          & pData$Age == 7] <- colors()[614]
pData$bg <- make.transparent(pData$col, alpha = 200)

pData$shp[pData$Mated == "mated"] <- 21
pData$shp[pData$Mated == "unmated"] <- 24

pData$Genotype <- relevel(pData$Genotype, ref = "N2")
pData$fileName <- rownames(pData)


## Get gene info (downloaded from wormbase):
tab <- read.table("Data/gene_info.txt",
                  sep = ";",
                  stringsAsFactors = FALSE)
fData <- data.frame(gene.id = c(rownames(count.matrix)[1:4], gsub("gene_id ", "", tab$V1)),
                    gene.name = c(rownames(count.matrix)[1:4], gsub(" gene_name ", "", tab$V2)),
                    source = c(rownames(count.matrix)[1:4], gsub(" gene_source ", "", tab$V3)),
                    biotype = c(rownames(count.matrix)[1:4], gsub(" gene_biotype ", "", tab$V4)),
                    stringsAsFactors = FALSE)
rownames(fData) <- fData$gene.name
identical(fData$gene.id, rownames(count.matrix))
## TRUE


## Create an expression set matrix to hold all the data and sample information:
eset <- ExpressionSet(assayData = count.matrix)
pData(eset) <- pData
fData(eset) <- fData

rownames(eset) <- fData(eset)$gene.name
colnames(eset) <- pData(eset)$Sample.Name

eset <- eset[, order(colnames(eset))]

## Remove the first 4 rows (unmapped reads etc):
eset <- eset[-c(1:4),]

## Filter out genes with low coverage across samples.
## Only keep genes that have a cpm (counts per million) value of at least 1 in at least 3 samples:
cov.per.sample <- colSums(exprs(eset))
norm.fact <- rep(cov.per.sample, each = nrow(eset))
cpm <- exprs(eset) / norm.fact * 10^6
ind.keep <- rowSums(cpm >= 1) >= 3
table(ind.keep)
## FALSE  TRUE
## 32601 14147
genes.keep <- names(ind.keep)[ind.keep]
eset <- eset[genes.keep, ]

## Gene type counts and eset dimensions:
table(fData(eset)$biotype)
##        antisense        lincRNA          miRNA          ncRNA         protein_coding     pseudogene 
##            3             33              1            60                     13743            286
##        rRNA         snoRNA          tRNA 
##         5            15               1 
dim(eset)
## Features  Samples
##     14147       32

## Coverage per sample:
quantile(colSums(exprs(eset)))
##       0%      25%      50%      75%     100%
##  8107127  9985455 10914898 12028495 13112002

## Number of genes per sample:
quantile(colSums(exprs(eset) != 0))
##      0%     25%     50%     75%    100%
##  13649.0 13778.5 13966.0 14049.5 14125.0

save(eset,
     file = "Data/eset_counts.RData")
