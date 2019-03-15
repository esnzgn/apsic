rm(list=ls())
source("waterfall_plot_methods.r")
source("common.r")


folder = ""
load(paste0(folder, "cancerData.RData"))


gene = "TP53"
waterfallForGene(cancerData, gene = gene, title=gene, rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="all", sig_alpha = NA)

waterfallForGene(cancerData, gene = gene, title=gene, rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_mut", sig_alpha = NA)

waterfallForGene(cancerData, gene = gene, title=gene, rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_wt", sig_alpha = NA)

waterfallForGene(cancerData, gene = gene, title=gene, rank=FALSE, legenedPos="bottomleft", 
                 cols=NULL, type="all", sig_alpha = NA)


selectedData = selectCelllines(cancerData, "Breast:Carcinoma")
gene = "TP53"
# copy number
waterfallForGene_CNA(selectedData, gene = gene, title = gene, rank=TRUE, legenedPos="bottomleft", cols=NULL, type="all", sig_alpha = NA) 

# mutation
waterfallForGene(selectedData, gene = gene, title = gene, rank=TRUE, legenedPos="bottomleft", cols=NULL, type="all", sig_alpha = NA) 


# gene names
gene_names = rownames(cancerData$viabilities)

# cancer names
cancer_types


# box plot
source("common.r")
library(edgeR)
library(ggplot2)


load("tcga/TCGA-BRCA.RData")
gene = "KRAS"
boxplot_gene(tcga_data, gene)


## cancer to TCGA project
cancerToTCGA_data = read.csv("CelllinesToTumorTypes.csv", stringsAsFactors = FALSE)
rownames(cancerToTCGA_data) = cancerToTCGA_data$Celllines

cancer_type = cancer_types[1]
cancerToTCGA_data[cancer_type, "mappedTCGAFile"]

cancer_type = cancer_types[10]
cancer_type
tmp = cancerToTCGA_data[cancer_type, "mappedTCGAFile"]
paste0("tcga/", tmp, ".RData")

load(paste0("tcga/", tmp, ".RData"))
gene = "KRAS"
boxplot_gene(tcga_data, gene)
