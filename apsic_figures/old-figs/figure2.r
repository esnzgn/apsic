rm(list=ls())
library(pheatmap)
library(reshape2)
library(ggplot2)
library(scales)

set.seed(1000)
driverGenesPanCancer <-  function() {
  cancer_type = "Pan_cancer"
  
  # missense
  dat_missense = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                                 "-missense-low.csv"), stringsAsFactors = FALSE, row.names = 1)
  n = sum(dat_missense$freq_mut > 10)
  gene_set3 = dat_missense[which(dat_missense$pvalue_mut < 1/n & dat_missense$pvalue_wt > 0.05), "gene"]
  
  # gene_set3 = dat_missense$gene[sample(nrow(dat_missense), 30)]
  
  # amplification
  dat_amplification = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                                      "-amplification-low.csv"), stringsAsFactors = FALSE, row.names = 1)
  n = sum(dat_amplification$freq_mut > 10)
  gene_set4 = dat_amplification[which(dat_amplification$pvalue_mut < 1/n & dat_amplification$pvalue_wt > 0.05), "gene"]
  # gene_set4 = gene_set3
  
  # truncating
  dat_truncating = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                                   "-truncating-high.csv"), stringsAsFactors = FALSE, row.names = 1)
  n = sum(dat_truncating$freq_mut > 10)
  gene_set5 = dat_truncating[which(dat_truncating$pvalue_mut < 1/n & dat_truncating$pvalue_wt > 0.05), "gene"]
  # gene_set5 = gene_set3
  
  gene_names = unique(c(gene_set3, gene_set4, gene_set5) )
  
  pvalues = cbind(dat_missense[gene_names, "pvalue_mut"], dat_truncating[gene_names, "pvalue_mut"],
                  dat_amplification[gene_names, "pvalue_mut"])
  pvalues[is.na(pvalues)] = 1
  colnames(pvalues) = c("missense", "truncating", "amplification")
  rownames(pvalues) = gene_names
  pvalues[order(pvalues[, 1]),]
}




prepareTCGADataMatrix <- function(direction) {
  pvalues = driverGenesPanCancer()
  folder = "../apsic_shiny/tcga/"
  fnames = list.files(folder)
  
  fname = fnames[1]
  gene_names = rownames(pvalues)
  for(fname in fnames) {
    fullName = paste0(folder, fname)

    load(fullName)
    print(paste(fname, ncol(tcga_data$tumor_count), ncol(tcga_data$normal_count)))
          
    if(direction == "low") {
      pvalues = cbind(pvalues, tcga_data$pvalues[gene_names, 1])  
    } else{
      pvalues = cbind(pvalues, tcga_data$pvalues[gene_names, 2])  
    }
  }
  
  pvalues[is.nan(pvalues)] = 1.0
  
  tmp = unname(sapply(fnames, function(x) { 
    x = substring(x, 6 )
    substring(x, 0, nchar(x)-6)} ))
  names(tmp) = tmp
  tmp["COADTCGA-READ"] = "COAD/READ"
  tmp["KIRCTCGA-KIRP"] = "KIRC/KIRP"
  
  colnames(pvalues)[4:ncol(pvalues)] = tmp
  pvalues
}

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

head(data_low)

data_low = prepareTCGADataMatrix("low")
plotHeatmap(data_low)

data_high = prepareTCGADataMatrix("high")

pdf("figures/fig2.pdf")
plotHeatmap(data_high)
dev.off()
