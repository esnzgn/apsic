rm(list=ls())

library(karyoploteR)


plotInChromosomeContext <- function(gene_names, geneAnnot) {


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
  
  
  # ids = setdiff(seq(1, n, 1), cosmic_indexes)
  # kpPlotMarkers(kp, data=genes[ids, ], labels=genes$hgnc_symbol[ids], r1 = 0.6, cex=0.75, text.orientation = "horizontal", data.panel = 1)
  
  kpPlotMarkers(kp, data=genes, labels=genes$gene_name, r1 = 0.6, cex=0.75, text.orientation = "horizontal", data.panel = 1)
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


# dat = read.csv("../apsic_shiny/apsic_pvalues/Breast_Carcinoma/Breast_Carcinoma-non-genetic-high.csv", stringsAsFactors = FALSE)
dat = read.csv("../apsic_shiny/apsic_pvalues/Breast_Carcinoma/Breast_Carcinoma-non-genetic-low.csv", stringsAsFactors = FALSE)
dat = read.csv("../apsic_shiny/apsic_pvalues/Breast_Carcinoma/Breast_Carcinoma-missense-low.csv", stringsAsFactors = FALSE)
dat = read.csv("../apsic_shiny/apsic_pvalues/Breast_Carcinoma/Breast_Carcinoma-amplification-low.csv", stringsAsFactors = FALSE)
dat = read.csv("../apsic_shiny/apsic_pvalues/Breast_Carcinoma/Breast_Carcinoma-truncating-high.csv", stringsAsFactors = FALSE)
n = nrow(dat)
# apsic_result = dat[which(dat$tumor_expressed_less < 0.05 ),]
apsic_result = dat[which(dat$pvalue_wt > 0.05 ),]
head(apsic_result)
gene_names = apsic_result$gene[1:40]

cancer_type = "Breast_Carcinoma"
gene_names = candidateGenes(cancer_type)

plotInChromosomeContext(unlist(gene_names), geneAnnot)
