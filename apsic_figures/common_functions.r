library(pheatmap)
library(reshape2)
library(ggplot2)
library(scales)
library(Rtsne)
library(factoextra)
library(metap)


allGenes <- function(cancer) {
  feature_type = "non-genetic-high"
  gene_names = c()
  # for(cancer in cancers) {
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer, "/", cancer, 
                        "-",feature_type, ".csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
  
  gene_names = rownames(dat)
}



selectImportantNongeneticGenes <- function(feature_type, cancers, nrGenesPerCancer, removeDuplicates) {
  gene_names = c()
  for(cancer in cancers) {
    dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer, "/", cancer, 
                          "-",feature_type, ".csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
    
    gene_names = c(gene_names, head(dat[order(as.numeric(dat$pvalue_wt), decreasing = F), "gene" ], nrGenesPerCancer) )
  }
  if(removeDuplicates){
    return(unique(gene_names))  
  }
  
  return(gene_names)
}

prepareNongeneticPvalueData <- function(feature_type, cancers) {
  
  # select genes
  gene_names = allGenes(cancers[2])
  pvalues = matrix(0, length(gene_names), length(cancers))
  rownames(pvalues) = gene_names
  colnames(pvalues) = cancers
  i = 1
  for(cancer in cancers) {
    dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer, "/", cancer,
                          "-",feature_type, ".csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
    pvalues[, i] = dat[gene_names,"pvalue_wt"]
    
    i = i+1
  }
  # pvalues = pvalues[order(pvalues[, "Pan_cancer"], decreasing = FALSE), ]
  # pvalues[is.na(pvalues)] = 1.0
  pvalues
}



plot_cluster=function(data, var_cluster, palette)  
{
  ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
    geom_point(size=0.25) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle("") +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal") + 
    scale_colour_brewer(palette = palette) 
}



plotHeatmap <- function(dat) {
  dat = dat[rev(rownames(dat)),]
  
  dat.m = melt(dat)
  names(dat.m) <- c("gene", "type", "pvalue")
  dat.m$pvalue[dat.m$pvalue > 20] = 20
  handle = ggplot(dat.m, aes(type, gene),font=3) +
    geom_tile(aes(fill = pvalue), color = "white") + 
    scale_fill_gradientn(name = expression(-log[10](p-value)), values=rescale(c(0, log10(0.05), 3, 5, 10, 20)), colours= c( "white",  "white",  "orange", "red1", "red2", "black")) +
    ylab("") +
    xlab("") +
    theme_bw()+
    theme(plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.y = element_text(size=7, face="italic"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
  
  
  print(handle)
  
}


