rm(list=ls())
library(pheatmap)
library(reshape2)
library(ggplot2)
library(scales)

plotHeatmap <- function(dat) {
  dat = dat[rev(rownames(dat)),]
  
  dat.m = melt(dat)
  names(dat.m) <- c("gene", "type", "pvalue")
  # dat.m$pvalue[dat.m$pvalue < 10^(-20)] = 10^(-20)
  handle = ggplot(dat.m, aes(type, gene),font=3) +
    geom_tile(aes(fill = pvalue), color = "white") + 
    # scale_fill_gradientn(values=rescale(c(10^(-20), 10^(-5), 0.001, 0.01, 0.3,  1)), colours=c("red3", "red2", "red1", "orange", "white",  "white") ) +
    scale_fill_gradientn(values=rescale(c(-20, -5, -3, -2, log10(0.05),  0)), colours=c("red3", "red2", "red1", "orange", "white",  "white") ) +
    ylab("") +
    xlab("") +
    theme_bw()+
    theme(plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "top")
  
  
  print(handle)

}

selectImportantGenesForFeature <- function(feature_type) {
  cancer_type = "Pan_cancer"
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                                 "-", feature_type, ".csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
  dat <- data.frame(dat)
  
  dat_sort= dat[order(as.numeric(dat$pvalue_mut)),]
  
  n = length(which(dat_sort$pvalue_mut <0.05))
  
  dat_sort[1:min(50, n) , "gene"]
}

preparePvalueData <- function(feature_type, cancers) {
  
  # select genes
  gene_names = selectImportantGenesForFeature(feature_type)
  # cancer = "Bladder_Carcinoma"
  pvalues = matrix(0, length(gene_names), length(cancers))
  rownames(pvalues) = gene_names
  colnames(pvalues) = cancers
  i = 1
  for(cancer in cancers) {
    dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer, "/", cancer, 
                          "-",feature_type, ".csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
    pvalues[, i] = dat[gene_names,"pvalue_mut"]
    
    i = i+1
  }
  pvalues = pvalues[order(pvalues[, "Pan_cancer"], decreasing = FALSE), ]
  # pvalues[is.na(pvalues)] = 1.0  
  pvalues
}

selectImportantNongeneticGenes <- function(feature_type, cancers) {
  gene_names = c()
  for(cancer in cancers) {
    dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer, "/", cancer, 
                          "-",feature_type, ".csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
    
    # print(cancer)
    # print(head(dat[order(as.numeric(dat$pvalue_wt), decreasing = F), "gene" ], 5) )
    
    gene_names = c(gene_names, head(dat[order(as.numeric(dat$pvalue_wt), decreasing = F), "gene" ], 5) )
  }
  unique(gene_names)
  # gene_names
}

prepareNongeneticPvalueData <- function(feature_type, cancers) {
  
  # select genes
  gene_names = selectImportantNongeneticGenes(feature_type, cancers)
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


cancers = list.files("../apsic_shiny/apsic_pvalues/")
cancers = c(cancers[which(cancers == "Pan_cancer")],cancers[-which(cancers == "Pan_cancer")])

feature_type = "missense-low"

pvalues = preparePvalueData("missense-low", cancers)
plotHeatmap(pvalues)

pvalues = preparePvalueData("amplification-low", cancers)
plotHeatmap(pvalues)

pvalues = preparePvalueData("truncating-high", cancers)
plotHeatmap(pvalues)



library(metap)
pvalues = prepareNongeneticPvalueData("non-genetic-high", cancers[-which(cancers == "Pan_cancer")])
combined_pvalues = apply(pvalues, 1, function(x) {
  res = sumlog(x[is.na(x) == FALSE])
  res$p
})

pvalues = pvalues[order(combined_pvalues), ]

pdf("non-genetic-high.pdf", 7, 15)
plotHeatmap(log10(pvalues) )
dev.off()



pvalues = prepareNongeneticPvalueData("non-genetic-low", cancers[-which(cancers == "Pan_cancer")])
combined_pvalues = apply(pvalues, 1, function(x) {
  res = sumlog(x[is.na(x) == FALSE])
  res$p
})

pvalues = pvalues[order(combined_pvalues), ]

pdf("non-genetic-low.pdf", 7, 17)
plotHeatmap(log10(pvalues) )
dev.off()
