rm(list=ls())
source("waterfall_plot_methods.r")
source("common.r")


folder = "ForShiny/"
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
gene = "CDK4"
waterfallForGene_CNA(selectedData, gene = gene, title = gene, rank=TRUE, legenedPos="bottomleft", cols=NULL, type="all", sig_alpha = NA) 

