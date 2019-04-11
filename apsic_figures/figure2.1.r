rm(list=ls())
library(pheatmap)
library(reshape2)
library(ggplot2)
library(scales)



plotHeatmap <- function(dat) {
  
  dat = dat[rev(rownames(dat)),]
  
  dat.m = melt(dat)
  names(dat.m) <- c("gene", "type", "pvalue")
  dat.m$pvalue[dat.m$pvalue < 10^(-20)] = 10^(-20)
  handle = ggplot(dat.m, aes(type, gene)) +
    geom_tile(aes(fill = pvalue), color = "white") + 
    scale_fill_gradientn(values=rescale(c(10^(-20), 10^(-5), 0.001,  1)), colours=c("red3", "red2", "red1",  "white") ) +
    ylab("") +
    xlab("") +
    theme_bw()+
    theme(plot.title = element_text(size=16),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") 
  print(handle)
}

cancer_type = "Pan_cancer"

dat_missense = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                               "-missense-low.csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
dat_missense <- data.frame(dat_missense)
typeof(dat_missense)

dat_missense_sort= dat_missense[order(as.numeric(dat_missense$pvalue_mut)),]

# head(dat_missense_sort)
# typeof(dat_missense)
# class(dat_missense)

gene_names = dat_missense_sort[1:30, "gene"]
# dat = read.csv("../apsic_shiny/apsic_pvalues/Bladder_Carcinoma/Bladder_Carcinoma-missense-low.csv", row.names = 1)
# # head(dat)
# 
# 
# mis_bladder <- dat[gene_names,"pvalue_wt"]


########## gene names & cancer

f_pvalue_vector <- function(gene_names, cancer){
  
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer, "/", cancer, 
                                 "-missense-low.csv"), stringsAsFactors = FALSE, row.names = 1,header = TRUE)
   dat[gene_names,"pvalue_wt"]
  
}
# ######

cancers = list.files("../apsic_shiny/apsic_pvalues/")

cancers = cancers[-which(cancers == "Pan_cancer")]

# cancer = "Bladder_Carcinoma"
pvalues = matrix(0, 30, length(cancers))
rownames(pvalues) = gene_names
colnames(pvalues) = cancers
i = 1
for(cancer in cancers) {
  pvalues[, i] = f_pvalue_vector(gene_names,cancer)  
   i = i+1
}
plotHeatmap(pvalues)





