rm(list=ls())

library(karyoploteR)


plotInChromosomeContext <- function(candidGenes, geneAnnot) {
  # TODO solve duplicated genes
  
  colors = c(rep("red", length(candidGenes$nongen_low)) , rep("blue", length(candidGenes$nongen_high)), 
             rep("yellow", length(candidGenes$gen_missense)), rep("green", length(candidGenes$gen_amplification)),
             rep("purple", length(candidGenes$gen_truncating)))
  
  
  gene_names = names(colors) = unname(unlist(candidGenes))
  
  
  # gene_names = unique(unlist(candidGenes) )
    # compute list of chromosomes which contain gene_names  
  genes_hg38_ = geneAnnot[which(geneAnnot$gene_name %in% gene_names), c("seqnames", "start", "end", "gene_name", "width")]
  genes_hg38 = genes_hg38_[order(genes_hg38_$gene_name, genes_hg38_$width), 1:4]
  
  rownames(genes_hg38) = NULL
  
  genes <- toGRanges(genes_hg38)
  seqlevelsStyle(genes) <- "UCSC"
  
  genes = genes[order(genes@seqnames, genes@ranges@start), ]
  
  char_element = as.character(substring(unique(as.character(genes@seqnames)), 4))
  
  allChromo = c(1:22, "X", "Y")
  u_chros = allChromo[which(allChromo %in% char_element)]
  

  kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes=paste0("chr",u_chros) )

  n = length(genes@seqnames)
  
  
  ids = seq(2, length(genes), by=2)
  
  tmp_genes = sapply(genes$gene_name[ids], toString)
  
  colors2 = colors[tmp_genes]
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 1, label.color=colors2)
  
  ids = seq(1, length(genes), by=2)
  colors2 = colors[tmp_genes]
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color=colors2)
  
  
  # 
  # if(length(candidGenes$nongen_low) > 0)
  # {
  #  
  #   ids = which(genes$gene_name %in% candidGenes$nongen_low)
  #   
  #   colors2 = colors[genes$gene_name[ids]]
  #   kpPlotMarkers(kp, data=genes[ids, ], labels=genes$gene_name[ids], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color=colors2)
  # }
  
# 
#   if(length(candidGenes$nongen_high) > 0)
#   {
#     ids = which(genes$gene_name %in% candidGenes$nongen_high)
#     kpPlotMarkers(kp, data=genes[ids, ], labels=genes$gene_name[ids], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 1, label.color="#d7191c")
#   }
#   
#   if(length(candidGenes$gen_missense) > 0)
#   {
#     ids = which(genes$gene_name %in% candidGenes$gen_missense)
#     kpPlotMarkers(kp, data=genes[ids, ], labels=genes$gene_name[ids], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color="#d7191c")
#   }
#   
#   if(length(candidGenes$gen_amplification) > 0)
#   {
#     ids = which(genes$gene_name %in% candidGenes$gen_amplification)
#     kpPlotMarkers(kp, data=genes[ids, ], labels=genes$gene_name[ids], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 1, label.color="#d7191c")
#   }
#   
#   if(length(candidGenes$gen_truncating) > 0)
#   {
#     ids = which(genes$gene_name %in% candidGenes$gen_truncating)
#     kpPlotMarkers(kp, data=genes[ids, ], labels=genes$gene_name[ids], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color="#2c7bb6")
#   }
  
  
  # evens = seq(2, length(genes), by=2)
  # kpPlotMarkers(kp, data=genes[evens, ], labels=genes$gene_name[evens], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2)
}




candidateGenes  <- function(cancer_type ) {
  
  # non-genetic high : potential tumor suppressor
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                 "-non-genetic-high.csv"), stringsAsFactors = FALSE)
  n = nrow(dat)
  gene_set1 = dat[which(dat$pvalue_wt < 1/n & dat$tumor_expressed_less <0.05), "gene"]

  # non-genetic low
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-non-genetic-low.csv"), stringsAsFactors = FALSE)
  gene_set2 = dat[which(dat$pvalue_wt < 1/n & dat$tumor_expressed_more <0.05), "gene"]
  
  # missense
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-missense-low.csv"), stringsAsFactors = FALSE)
  n = sum(dat$freq_mut > 1)
  gene_set3 = dat[which(dat$pvalue_mut < 1/n & dat$pvalue_wt > 0.05), "gene"]

  # amplification
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-amplification-low.csv"), stringsAsFactors = FALSE)
  n = sum(dat$freq_mut > 1)
  gene_set4 = dat[which(dat$pvalue_mut < 1/n & dat$pvalue_wt > 0.05), "gene"]
  
  
  # truncating
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-truncating-high.csv"), stringsAsFactors = FALSE)
  n = sum(dat$freq_mut > 1)
  gene_set5 = dat[which(dat$pvalue_mut < 1/n & dat$pvalue_wt > 0.05), "gene"]
  
  
  list(nongen_high=gene_set1, nongen_low=gene_set2, gen_missense=gene_set3, 
       gen_amplification=gene_set4, gen_truncating=gene_set5)
}

#############################

geneAnnot = read.table("Homo_sapiens.GRCh38.p10_genes.csv", sep='\t', stringsAsFactors = FALSE)

cancer_type = "Breast_Carcinoma"
# cancer_type = "Liver_HCC"
candidGenes = candidateGenes(cancer_type)

# gene_names = unname(unlist(candidGenes))
# typeOfGenes

plotInChromosomeContext(candidGenes, geneAnnot)
# 
# load("../apsic_shiny/cancerData.RData")
# a = 0:4
# for (i in a){
#   i = i + 1
#   print(cancer_type)
#   pdf(paste0("KP", cancer_type)) <- plotInChromosomeContext(candidGenes, geneAnnot)
#   
# 
# }
# cancer_type
