# Objective: Perform Pathway analysis on DEGs (FDR < 0.05) using Fischer's exact test
# Use DEG table (with fold change and FDR-values) created by DEseq2  


rm(list=ls())

# Set working directory
setwd("~/Dropbox/Code-For-Reviewers/")

genes <- read.delim("Results/DEseq/Groups/N3MvF3M.txt", 
                      stringsAsFactors = FALSE,
                      row.names = 1)

# Remove N/A rows
genes <- na.omit(genes)

# Upload GO terms (note: the script will call these kegg but they are GO terms)
load(file = "Data/Cele_Go_terms.rda")

########################################################################
# Extract genes that are DOWN in the comparison (fem-1 enriched)
deg.down <- genes[which(genes$padj < 0.05 & genes$log2FoldChange < 0), ] #864 fem-1 enriched

# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.down <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.down) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.down) <- colnames(gene.set.data)


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
  deg_kegg_genes <- nrow(genes[which(rownames(deg.down) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(deg.down[which(rownames(deg.down) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(deg.down) - deg_kegg_genes 
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg.down)) - non_deg_kegg_genes 
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data.down$p.values[i] <- ftest$p.value
  my.data.down$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.down$FDR.BH <- p.adjust(my.data.down$p.values, "BH")
# Order by FDR adjusted p-values
my.data.downOrdered <- my.data.down[order(my.data.down$FDR.BH),]
head(my.data.downOrdered, n=10)

# Write a table with the final results 
write.table(my.data.downOrdered, 
            col.names = NA,
            quote = FALSE,
            sep = "\t",
            file="Results/DEseq/Groups/GO/N3MvF3M-fem-1UP.txt")

#####################################################################################################
# Extract genes that are UP in N2
deg.up <- genes[which(genes$padj < 0.05 & genes$log2FoldChange > 0),] #1025 UP in N3M
 
# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.up <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.up) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.up) <- colnames(gene.set.data)
 

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
 deg_kegg_genes <- nrow(genes[which(rownames(deg.up) %in% gene.set.data[,i]),])
 deg_kegg_genes.names <- as.data.frame(row.names(deg.up[which(rownames(deg.up) %in% gene.set.data[,i]),]))
 deg_non_kegg_genes <- nrow(deg.up) - deg_kegg_genes 
 non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
 non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg.up)) - non_deg_kegg_genes 
 ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
 my.data.up$p.values[i] <- ftest$p.value
 my.data.up$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.up$FDR.BH <- p.adjust(my.data.up$p.values, "BH")
# Order by FDR adjusted p-values
my.data.upOrdered <- my.data.up[order(my.data.up$FDR.BH),]
head(my.data.upOrdered, n=10)
 
# Write a table with the final results 
write.table(my.data.upOrdered, 
            col.names = NA,
            quote = FALSE,
            sep = "\t",
            file="Results/DEseq/Groups/GO/N3MvF3M-N2UP.txt")
