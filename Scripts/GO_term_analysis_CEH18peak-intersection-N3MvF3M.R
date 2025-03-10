# Objective: Perform Pathway analysis on DEGs (FDR < 0.05) using Fischer's exact test
# Use DEG table (with fold change and FDR-values) created by DEseq2  
# Will use KEGG pathways from Param


rm(list=ls())

# Set working directory
setwd("~/Dropbox/Code-For-Reviewers/")

genes <- read.delim("Results/DEseq/Groups/N3MvF3M.txt", 
                      stringsAsFactors = FALSE,
                      row.names = 1)

CEH18.all <- readLines(con = "Data/CEH-18-peak-associations.txt")

# Remove N/A rows
genes <- na.omit(genes)

# Upload GO terms (note: the script will call these kegg but they are GO terms)
load(file = "Data/Cele_Go_terms.rda")

########################################################################
# Extract genes that are expressed more highly in N2 and that are associated with a CEH-18 peak
N2 <- genes[which(genes$padj < 0.05 & genes$log2FoldChange > 0), ] 
CEH18.N <- N2[which(rownames(N2) %in% CEH18.all), ] 

# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data) <- c("p.values", "FDR.BH", "genes")
rownames(my.data) <- colnames(gene.set.data)

## In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in my gene list (background) are found in the KEGG pathway (kegg_genes)
# 2) how many DEGs are found in the KEGG pathway (deg_kegg_genes)
# 3) how many DEGs are not found in the KEGG pathway (deg_non_kegg_genes)
# 4) how many genes in the KEGG pathway are not DEG (non_deg_kegg_genes)
# 5) how many genes are neither DEGs nor found in the KEGG pathway (non_deg_non_kegg_genes)
# 
# Note: deg_kegg_genes + non_deg_kegg_genes = kegg_genes

# Use the numbers to perform fischer's t-test, and save the p-value in the data frame 
# Finally correct for FDR using benjimini_hochberg

for(i in 1:ncol(gene.set.data)){
  kegg_genes <- nrow(genes[which(rownames(genes) %in% gene.set.data[,i]),])
  deg_kegg_genes <- nrow(genes[which(rownames(CEH18.N) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(CEH18.N[which(rownames(CEH18.N) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(CEH18.N) - deg_kegg_genes 
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(CEH18.N)) - non_deg_kegg_genes 
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data$p.values[i] <- ftest$p.value
  my.data$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data$FDR.BH <- p.adjust(my.data$p.values, "BH")
# Order by FDR adjusted p-values
my.dataOrdered <- my.data[order(my.data$FDR.BH),]
head(my.dataOrdered, n=10)

# Write a table with the final results 
write.table(my.dataOrdered, 
            col.names = NA,
            quote = FALSE,
            sep = "\t",
            file="Results/DEseq/Groups/GO/N3MvF3M-N3Mup-withCEH18peak.txt")



########################################################################
# Extract genes that are expressed more highly in fem-1 and that are associated with a CEH-18 peak
fem <- genes[which(genes$padj < 0.05 & genes$log2FoldChange < 0), ]
CEH18.F <- fem[which(rownames(fem) %in% CEH18.all), ]

# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data) <- c("p.values", "FDR.BH", "genes")
rownames(my.data) <- colnames(gene.set.data)

## In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in my gene list (background) are found in the KEGG pathway (kegg_genes)
# 2) how many DEGs are found in the KEGG pathway (deg_kegg_genes)
# 3) how many DEGs are not found in the KEGG pathway (deg_non_kegg_genes)
# 4) how many genes in the KEGG pathway are not DEG (non_deg_kegg_genes)
# 5) how many genes are neither DEGs nor found in the KEGG pathway (non_deg_non_kegg_genes)
# 
# Note: deg_kegg_genes + non_deg_kegg_genes = kegg_genes

# Use the numbers to perform fischer's t-test, and save the p-value in the data frame 
# Finally correct for FDR using benjimini_hochberg

for(i in 1:ncol(gene.set.data)){
  kegg_genes <- nrow(genes[which(rownames(genes) %in% gene.set.data[,i]),])
  deg_kegg_genes <- nrow(genes[which(rownames(CEH18.F) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(CEH18.F[which(rownames(CEH18.F) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(CEH18.F) - deg_kegg_genes 
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(CEH18.F)) - non_deg_kegg_genes 
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data$p.values[i] <- ftest$p.value
  my.data$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data$FDR.BH <- p.adjust(my.data$p.values, "BH")
# Order by FDR adjusted p-values
my.dataOrdered <- my.data[order(my.data$FDR.BH),]
head(my.dataOrdered, n=10)

# Write a table with the final results 
write.table(my.dataOrdered, 
            col.names = NA,
            quote = FALSE,
            sep = "\t",
            file="Results/DEseq/Groups/GO/N3MvF3M-F3Mup-withCEH18peak.txt")

