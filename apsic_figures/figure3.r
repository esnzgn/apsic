rm(list=ls())
library(karyoploteR)

candidGenes_redundency_scanner <- function(candidGenes){
  
  # nongene_high as template
  no_ngh = length(candidGenes$nongen_high)
  if (no_ngh > 0){
    for (i in no_ngh){
      no_ngl = length(candidGenes$nongen_low)
      for (z in no_ngl){
        if (candidGenes$nongen_high[i] == candidGenes$nongen_low[z]){
          candidGenes$nongen_high[i] = paste0(candidGenes$nongen_high[i], "*")}}
      z = 1
      no_gm = length(candidGenes$gen_missense)
      for (z in no_gm){
        if (candidGenes$nongen_high[i] == candidGenes$gen_missense[z]){
          candidGenes$nongen_high[i] = paste0(candidGenes$nongen_high[i], "@")}}
      z = 1
      no_ga = length(candidGenes$gen_amplification)
      for (z in no_ga){
        if (candidGenes$nongen_high[i] == candidGenes$gen_amplification[z]){
          candidGenes$nongen_high[i] = paste0(candidGenes$nongen_high[i], "$")}}
      z = 1
      no_gt = length(candidGenes$gen_truncating)
      for (z in no_gt){
        if (candidGenes$nongen_high[i] == candidGenes$gen_truncating[z]){
          candidGenes$nongen_high[i] = paste0(candidGenes$nongen_high[i], "&")}}
      
    }
  }
  # nongene_low as template
  no_ngl = length(candidGenes$nongen_low)
  if (no_ngl > 0){
    for (i in no_ngl){
      no_ngh = length(candidGenes$nongen_high)
      z = 1
      for (z in no_ngh){
        if (candidGenes$nongen_low[i] == candidGenes$nongen_high[z]){
          candidGenes$nongen_low[i] = paste0(candidGenes$nongen_low[i], "#")}}
      z = 1
      no_gm = length(candidGenes$gen_missense)
      for (z in no_gm){
        if (candidGenes$nongen_high[i] == candidGenes$gen_missense[z]){
          candidGenes$nongen_low[i] = paste0(candidGenes$nongen_low[i], "@")}}
      z = 1
      no_ga = length(candidGenes$gen_amplification)
      for (z in no_ga){
        if (candidGenes$nongen_low[i] == candidGenes$gen_amplification[z]){
          candidGenes$nongen_low[i] = paste0(candidGenes$nongen_low[i], "$")}}
      z = 1
      no_gt = length(candidGenes$gen_truncating)
      for (z in no_gt){
        if (candidGenes$nongen_low[i] == candidGenes$gen_truncating[z]){
          candidGenes$nongen_low[i] = paste0(candidGenes$nongen_low[i], "&")}}
      
    }
  }
  # gen_missense as template
  no_gm = length(candidGenes$gen_missense)
  if (no_gm > 0){
    for (i in no_gm){
      no_ngh = length(candidGenes$nongen_high)
      z = 1
      for (z in no_ngh){
        if (candidGenes$gen_missense[i] == candidGenes$nongen_high[z]){
          candidGenes$gen_missense[i] = paste0(candidGenes$gen_missense[i], "#")}}
      z = 1
      no_ngl = length(candidGenes$nongen_low)
      for (z in no_ngl){
        if (candidGenes$gen_missense[i] == candidGenes$nongen_low[z]){
          candidGenes$gen_missense[i] = paste0(candidGenes$gen_missense[i], "*")}}
      z = 1
      no_ga = length(candidGenes$gen_amplification)
      for (z in no_ga){
        if (candidGenes$gen_missense[i] == candidGenes$gen_amplification[z]){
          candidGenes$gen_missense[i] = paste0(candidGenes$gen_missense[i], "$")}}
      z = 1
      no_gt = length(candidGenes$gen_truncating)
      for (z in no_gt){
        if (candidGenes$gen_missense[i] == candidGenes$gen_truncating[z]){
          candidGenes$gen_missense[i] = paste0(candidGenes$gen_missense[i], "&")}}
      
    }
  }
  # gen_amplification as template
  no_ga = length(candidGenes$gen_amplification)
  if (no_ga > 0){
    for (i in no_ga){
      no_ngh = length(candidGenes$nongen_high)
      z = 1
      for (z in no_ngh){
        if (candidGenes$gen_amplification[i] == candidGenes$nongen_high[z]){
          candidGenes$gen_amplification[i] = paste0(candidGenes$gen_amplification[i], "#")}}
      z = 1
      no_ngl = length(candidGenes$nongen_low)
      for (z in no_ngl){
        if (candidGenes$gen_amplification[i] == candidGenes$nongen_low[z]){
          candidGenes$gen_amplification[i] = paste0(candidGenes$gen_amplification[i], "*")}}
      z = 1
      no_gm = length(candidGenes$gen_missense)
      for (z in no_ngl){
        if (candidGenes$gen_amplification[i] == candidGenes$gen_missense[z]){
          candidGenes$gen_amplification[i] = paste0(candidGenes$gen_amplification[i], "@")}}
      z = 1
      no_gt = length(candidGenes$gen_truncating)
      for (z in no_gt){
        if (candidGenes$gen_amplification[i] == candidGenes$gen_truncating[z]){
          candidGenes$gen_amplification[i] = paste0(candidGenes$gen_amplification[i], "&")}}
      
    }
  }
  # gen_truncating as template
  no_gt = length(candidGenes$gen_truncating)
  if (no_gt > 0){
    for (i in no_gt){
      no_ngh = length(candidGenes$nongen_high)
      z = 1
      for (z in no_ngh){
        if (candidGenes$gen_truncating[i] == candidGenes$nongen_high[z]){
          candidGenes$gen_truncating[i] = paste0(candidGenes$gen_truncating[i], "#")}}
      z = 1
      no_ngl = length(candidGenes$nongen_low)
      for (z in no_ngl){
        if (candidGenes$gen_truncating[i] == candidGenes$nongen_low[z]){
          candidGenes$gen_truncating[i] = paste0(candidGenes$gen_truncating[i], "*")}}
      z = 1
      no_gm = length(candidGenes$gen_missense)
      for (z in no_ngl){
        if (candidGenes$gen_truncating[i] == candidGenes$gen_missense[z]){
          candidGenes$gen_truncating[i] = paste0(candidGenes$gen_truncating[i], "@")}}
      z = 1
      no_ga = length(candidGenes$gen_amplification)
      for (z in no_ga){
        if (candidGenes$gen_truncating[i] == candidGenes$gen_amplification[z]){
          candidGenes$gen_truncating[i] = paste0(candidGenes$gen_truncating[i], "&")}}
      
    }
  } 
}


plotInChromosomeContext <- function(candidGenes, geneAnnot, plot_title) {
  # TODO solve duplicated genes

  colors = c(rep("red", length(candidGenes$nongen_low)) , rep("blue", length(candidGenes$nongen_high)), 
             rep("orange", length(candidGenes$gen_missense)), rep("green", length(candidGenes$gen_amplification)),
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
  

  kp <- plotKaryotype(genome="hg38", plot.type = 2, chromosomes=paste0("chr",u_chros), main= plot_title )

  n = length(genes@seqnames)
  
  
  ids = seq(2, length(genes), by=2)
  tmp_genes = sapply(genes$gene_name[ids], toString)
  
  colors2 = colors[tmp_genes]
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 1, label.color=colors2, font=3)
  
  ids = seq(1, length(genes), by=2)
  tmp_genes = sapply(genes$gene_name[ids], toString)
  colors2 = colors[tmp_genes]
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color=colors2, font=3)
  legend("bottomright", legend=c("Oncogene mutation", "Oncogene amplification", "Tumor suppressor mutation",
                                 "Non-genetic oncogene", "Non-genetic tumor suppressor"), 
         col=c("orange", "green", "purple", "blue", "red"), cex=1, lty=1, lwd=2)
}




candidateGenes  <- function(cancer_type ) {
  
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-missense-low.csv"), stringsAsFactors = FALSE)
  
  nrOfCelllines = dat[1, ]$freq_mut + dat[1,]$freq_wt
  minNrOfCelllines=ceiling(nrOfCelllines/100)
  
  gene_set1 = gene_set2 = ""
  # non-genetic high : potential tumor suppressor
  fname = paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                 "-non-genetic-high.csv")
  if(file.exists(fname)) {
    dat = read.csv(fname, stringsAsFactors = FALSE)
    n = sum(dat$freq_wt > minNrOfCelllines)
    gene_set1 = dat[which(dat$pvalue_wt < 1/n & dat$tumor_expressed_less <0.05), "gene"]
  }

  # non-genetic low
  fname = paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                 "-non-genetic-low.csv")
  if(file.exists(fname)) {
    dat = read.csv( fname, stringsAsFactors = FALSE)
    n = sum(dat$freq_wt > minNrOfCelllines)
    gene_set2 = dat[which(dat$pvalue_wt < 1/n & dat$tumor_expressed_more <0.05), "gene"]
  }
  
  # missense
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-missense-low.csv"), stringsAsFactors = FALSE)
  n = sum(dat$freq_mut > minNrOfCelllines)
  gene_set3 = dat[which(dat$pvalue_mut < 1/n & dat$pvalue_wt > 0.05), "gene"]

  # amplification
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-amplification-low.csv"), stringsAsFactors = FALSE)
  n = sum(dat$freq_mut > minNrOfCelllines)
  gene_set4 = dat[which(dat$pvalue_mut < 1/n & dat$pvalue_wt > 0.05), "gene"]
  
  
  # truncating
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-truncating-high.csv"), stringsAsFactors = FALSE)
  n = sum(dat$freq_mut > minNrOfCelllines)
  gene_set5 = dat[which(dat$pvalue_mut < 1/n & dat$pvalue_wt > 0.05), "gene"]
  
  
  list(nongen_high=gene_set1, nongen_low=gene_set2, gen_missense=gene_set3, 
       gen_amplification=gene_set4, gen_truncating=gene_set5)
}

#############################

geneAnnot = read.table("Homo_sapiens.GRCh38.p10_genes.csv", sep='\t', stringsAsFactors = FALSE)

cancer_type = "Breast_Carcinoma"
# cancer_type = "Liver_HCC"
candidGenes = candidateGenes(cancer_type)
candidGenes_redundency_scanner(candidGenes)
plotInChromosomeContext(candidGenes, geneAnnot, cancer_type)


output_folder =  "figures/fig3/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
file_names = list.files("../apsic_shiny/apsic_pvalues/")

for(fname in file_names)
{
  candidGenes = candidateGenes(fname)
  pdf(paste0(output_folder, fname, ".pdf"), 10, 15)
  plotInChromosomeContext(candidGenes, geneAnnot, fname)
  dev.off()
}



pdf(paste0(output_folder, "All_cancers", ".pdf"), 10, 15,  onefile=TRUE)
for(fname in file_names) {
  candidGenes = candidateGenes(fname)
  plotInChromosomeContext(candidGenes, geneAnnot, fname)
}

dev.off()
