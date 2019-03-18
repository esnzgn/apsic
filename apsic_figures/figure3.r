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
  tmp_genes = sapply(genes$gene_name[ids], toString)
  colors2 = colors[tmp_genes]
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color=colors2)
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

# cancer_type = "Breast_Carcinoma"
# # cancer_type = "Liver_HCC"
# candidGenes = candidateGenes(cancer_type)
# 

output_folder =  "figures/fig3/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
file_names = list.files("../apsic_shiny/apsic_pvalues/")

for(fname in file_names)
{
  candidGenes = candidateGenes(fname)
  pdf(paste0(output_folder, fname, ".pdf"), 14, 8)
  plotInChromosomeContext(candidGenes, geneAnnot)
  dev.off()
}
