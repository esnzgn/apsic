rm(list=ls())
library(pvclust)
source("scripts/common_functions.r")

plotClustersOfCancers <- function(pvalues, fig_folder, feature_name, nrGenes=500, nboot=1000) {
  ################# clustering cancer subtypes according to the p-values ###############################
  ######################################################################################################
  ### first select 500 genes with highly variable p-values among subtypes
  sorted_genes = sort(apply(pvalues, 1, var), decreasing = TRUE)
  genes = names(sorted_genes)[1:nrGenes]
  dat = pvalues[genes, ]
  
  # select NA values to 1, if any. (In our case, there is no NA p-value so the following command has no impact)
  dat[which(is.na(dat))] = 1.0
  
  # clustering with 10000 bootstraps 
  result = pvclust(dat, method.dist="cor", method.hclust="ward.D2", nboot=nboot)
  
  pdf(paste0(fig_folder, "pvclust_90_", feature_name, ".pdf"), width=10, height=8)
  plot(result, hang=-1, main='', print.num=FALSE, xlab='', sub='')
  pvrect(result, alpha=0.90)
  dev.off()
  
  as.dendrogram(result$hclust)
  pdf(paste0(fig_folder, "dendrogram_horiz_90_", feature_name, ".pdf"), width=3, height=20) 
  plot(rev(as.dendrogram(result$hclust)), type = "rectangle", ylab = "Height", horiz = TRUE, leaflab="none")
  dev.off()

  pdf(paste0(fig_folder, "pvclust_85_", feature_name, ".pdf"), width=10, height=8)
  plot(result, hang=-1, main='', print.num=FALSE, xlab='', sub='')
  pvrect(result,  alpha=0.85)
  dev.off()
  
  pdf(paste0(fig_folder, "dendrogram_horiz_85_", feature_name, ".pdf"), width=3, height=20) 
  plot(rev(as.dendrogram(result$hclust)), type = "rectangle", ylab = "Height", horiz = TRUE, leaflab="none")  
  dev.off()

  pdf(paste0(fig_folder, "pvclust_80_", feature_name, ".pdf"), width=10, height=8)
  plot(result, main='', hang=-1, print.num=FALSE, xlab='', sub='')
  pvrect(result, alpha=0.80)
  dev.off()
  
  pdf(paste0(fig_folder, "dendrogram_horiz_80_", feature_name, ".pdf"), width=3, height=20) 
  plot(rev(as.dendrogram(result$hclust)), type = "rectangle", ylab = "Height", horiz = TRUE, leaflab="none")
  dev.off()

  pdf(paste0(fig_folder, "dendrogram_horiz_80_", feature_name, "_with_labs.pdf"), width=3, height=20) 
  plot(rev(as.dendrogram(result$hclust)), type = "rectangle", ylab = "Height", horiz = TRUE)
  dev.off()
  
  # now we reorder subtypes according to the hierarchical clustering
  clusters = hclust(get_dist(t(dat), method = "pearson"), method= "ward.D2")
  plot(clusters)
  cancers = clusters$labels[clusters$order]
  save(cancers, file=paste0(fig_folder2, "ordered_cancers_",feature_name, ".RData"))
}

fig_folder2 = "figures/figS1/"
dir.create(fig_folder2, recursive = T, showWarnings = FALSE)

##### we first load all p-value files
cancers = list.files("../apsic_shiny/apsic_pvalues/")
cancers = cancers[-which(cancers == "Pan_cancer")]

set.seed(10)
pvalues = prepareNongeneticPvalueData("non-genetic-tumor_suppressor", cancers)
plotClustersOfCancers(pvalues, fig_folder2, "tumor_suppressor", nboot=10000)


set.seed(10)
pvalues = prepareNongeneticPvalueData("non-genetic-oncogene", cancers)
plotClustersOfCancers(pvalues, fig_folder2, "oncongenes", nboot=10000)
