rm(list=ls())

library(karyoploteR)


plotInChromosomeContext <- function(gene_names, geneAnnot) {

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


  odds = seq(1, length(genes), by=2)
  kpPlotMarkers(kp, data=genes[odds, ], labels=genes$gene_name[odds], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 1, label.color="#7b3294")
  
  evens = seq(2, length(genes), by=2)
  kpPlotMarkers(kp, data=genes[evens, ], labels=genes$gene_name[evens], r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2)
}


geneAnnot = read.table("Homo_sapiens.GRCh38.p10_genes.csv", sep='\t')


dat = read.csv("../apsic_shiny/apsic_pvalues/Pan_cancer/Pan_cancer-missense-low.csv", stringsAsFactors = FALSE)
n = nrow(dat)
apsic_result = dat[which(dat$pvalue_wt > 0.05 & dat$pvalue_mut < 1/n),]
gene_names = apsic_result$gene

plotInChromosomeContext(gene_names, geneAnnot)


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

cancer_type = "Breast_Carcinoma"
# cancer_type = "Liver_HCC"
gene_names = candidateGenes(cancer_type)

plotInChromosomeContext(unique(unlist(gene_names) ), geneAnnot)
