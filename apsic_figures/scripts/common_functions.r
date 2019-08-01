library(pheatmap)
library(reshape2)
library(ggplot2)
library(scales)
library(Rtsne)
library(factoextra)
library(metap)


allGenes <- function(cancer) {
  feature_type = "non-genetic-tumor_suppressor"
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
    
    topGeneRows = head(dat[order(as.numeric(dat$pvalue), decreasing = F), ], nrGenesPerCancer)     
    gene_names = c(gene_names, rownames(topGeneRows))
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
    pvalues[, i] = dat[gene_names, "pvalue"]
    
    i = i+1
  }
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
#dat = -log10(pvalues)
plotHeatmap <- function(dat) {
  dat = dat[rev(rownames(dat)),]
  orgGeneNames = rownames(dat)
  dat = data.frame(dat)
  geneNames = rownames(dat) 
  dat$gene = geneNames
  names(orgGeneNames) = geneNames
  
  group = rev(rep(1:(nrow(dat)/5), each=5))
  names(group) = geneNames
  
  dat.m = melt(dat, id.vars = "gene")
  names(dat.m) <- c("gene", "type", "pvalue")
  
  mylevels <- geneNames
  dat.m$gene <- factor(dat.m$gene,levels=mylevels)
  
  dat.m$pvalue[dat.m$pvalue > 20] = 20
  
  dat.m$group = group[dat.m$gene]
  orgGenes  = orgGeneNames[dat.m$gene]
  
  handle = ggplot(dat.m, aes(type, gene),font=3) +
    facet_grid(group~., scales = "free_y") +
    geom_tile(aes(fill = pvalue), color = "gray") + 
    scale_fill_gradientn(name = expression(-log[10](p-value)), values=rescale(c(0, log10(0.05), 3, 5, 10, 20)), colours= c( "white",  "white",  "orange", "red1", "red2", "black")) +
    scale_y_discrete(labels=orgGenes) +
    ylab("") +
    xlab("") +
    theme_bw()+
    theme(plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.y = element_text(size=6, face="italic"),
          axis.text.x = element_text(size=6, angle = 45, hjust = 1),
          legend.position = "right")
  
  print(handle)
  
}