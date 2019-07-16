These scripts were used to analyze the RNA-seq data for the manuscript "Self-sperm promotes resistance to male-induced demise"

A description of the naming scheme for the data:
For each sample, there is a three position code that describes it. The first position is the genotype (N is WT N2 hermaphrodites, F is fem-1 feminized individuals). The second position is the age when the exposure to males occurred (day 3 or 7 of life). Note that the worms were collected for RNA-seq on the following day. The final position indicates their mating-status (U is for the individuals that never encountered males and M is for the individuals that successfully mated and received male sperm). For example: F3M are young feminized individuals that received male sperm following a brief interaction with males.

An overview of the data flow:
1. Import.R
	This script takes the read count values, sample information, and gene information to organize the data for DEseq2. Import.R filters the data so that genes with low coverage are removed. This script also calculates basic information about the datasets (coverage, gene type counts etc.). 

2. Normalize.R
	This script uses DEseq2 to calculate differential expression and generate a data dispersion plot for data quality testing. Normalize.R normalizes the datasets using variance stabilization transformation. This normalized data is used to generate PCA plots. Note that the script generates 3 sets of VST normalized data--the full dataset, only the RNA-seq data from the hermpahrodites that never encountered males, and only the RNA-seq data of the mated hermaphrodites that received male sperm.

3. Descriptive-*.R
	This script performs principal component analysis and generates PCA plots and a dendrogram of the data. There are three versions of this script for creating a PCA of all the data (Figure 2--figure supplement 1B) and for the two PCAs displayed in Fig. 2E, F. The colors shown in the manuscript differ slightly from the plots made here because the coloring was modified in Illustrator.

4. Model_with_Groups.R
	This script uses DEseq2 to calculate differential expression. The values calculated here are used in GO_term_analysis_*.R and are included in a supplemental table.

5. GO_term_analysis_*.R
	This set of scripts calculates GO term enrichment for the gene sets of interest. These scripts use the differential expression calculations from Model_with_Groups.R to generate the data displayed in Figure 2 and the supplemental data and the CEH-18 peak associated genes for Figure 4 and supplemental data. The bar plots were generated with Prism 7 using -log10 transformed data from this script. 

6. Heatmaps_Pathways_Read_Counts_*.R
	In this script, the appropriate sets of samples are subset and normalized by variance stabilization transformation. These normalized read counts are used to generate heatmaps that are displayed in Fig. 2, G and H. This script uses the GO gene lists to generate the subset of genes that are used for each heatmap. In the figures, several distinct heatmaps are combined displayed in each figure (separated by white space between the different gene sets). Note that some genes may appear more than once in a figure panel because they are members of more than one GO group.

(7. make.transparent.R is only used for aesthetics)

8. "pre_processing_DEGs.R" 
	Reads in RNA-seq DEG list and output those with an FDR<0.05 that are upregulated per condition
	
9. "get_DEG_promoters.R" 
	Uses the DEG lists from "pre_processing_DEGs.R" and outputs a bed file for each containing the associated promoter (definied as +/- 300bp around the TSS).

10. "intersect_CEH18_targets.R" 
	Outputs a list of CEH-18 ChIP-seq peaks that are associated with different lists of DEGs as well as outputs a unique gene list containing the DEGs that have nearby CEH-18 ChIP-seq peaks


Package and R Versions:
R version 3.2.4
DEseq2 version 1.10.1
Biobase version 2.30.0

Data sets used:
C. elegans genome: WBcel235
C. elegant gene ontology associations: validated on 04/15/2016 
Gene Ontology: releases/2016-06-02
GO from: www.geneontology.org
