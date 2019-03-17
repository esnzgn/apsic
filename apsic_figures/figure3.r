rm(list=ls())

library(karyoploteR)


plotInChromosomeContext <- function(gene_names, geneAnnot) {


  genes_hg38_ = geneAnnot[which(geneAnnot$gene_name %in% gene_names), c("seqnames", "start", "end", "gene_name", "width")]
  genes_hg38 = genes_hg38_[order(genes_hg38_$gene_name, genes_hg38_$width), 1:4]
  
  rownames(genes_hg38) = NULL
  
  genes <- toGRanges(genes_hg38)
  seqlevelsStyle(genes) <- "UCSC"
  
  genes = genes[order(genes@seqnames, genes@ranges@start), ]
  
  u_chros = unique(as.character(genes@seqnames))
  u_chros = u_chros[order(as.numeric(substring(u_chros, 4)))]
  
  
  kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes=u_chros)

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


dat = read.csv("../apsic_shiny/apsic_pvalues/Breast_Carcinoma/Breast_Carcinoma-non-genetic-high.csv", stringsAsFactors = FALSE)
n = nrow(dat)
apsic_result = dat[which(dat$pvalue_wt < 1/n & dat$tumor_expressed_less < 0.05 ),]
gene_names = apsic_result$gene

plotInChromosomeContext(gene_names, geneAnnot)
