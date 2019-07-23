rm(list=ls())
library(karyoploteR)

annotateGene <- function(candidGenes){
  
  ##################################################################################################
  ##################################### redundency Matrix scan #####################################
  ##################################################################################################
  
  genes = unique(unname(unlist(candidGenes)))
  
  x <- matrix(0, nrow = length(genes), ncol = 7)
  rownames(x) = genes
  colnames(x) <- c("nongen_high","nongen_low","gen_missense","gen_amplification","gen_truncating","color","label")
  x = data.frame(x)
  
  gene = 1
  for (gene in genes){
    if (gene %in% candidGenes$nongen_high){
      x[gene, "nongen_high"] = 1
    }
    if (gene %in% candidGenes$nongen_low){
      x[gene, "nongen_low"] = 1
    }
    if (gene %in% candidGenes$gen_missense){
      x[gene,"gen_missense"] = 1
    }
    if (gene %in% candidGenes$gen_amplification){
      x[gene,"gen_amplification"] = 1
    }
    if (gene %in% candidGenes$gen_truncating){
      x[gene,"gen_truncating"] = 1
    }
  }
  
  ############################################################################
  ########################### Color identification ###########################
  ############################################################################
  gene = 1
  for (gene in genes){
    x[gene,"label"] = rownames(x[gene,])
    if (sum(x[gene, "nongen_high"]+x[gene, "nongen_low"]+x[gene,"gen_missense"]+x[gene,"gen_amplification"]+x[gene,"gen_truncating"]) > 1){
      x[gene,"color"] = "black"
      # ,x[gene,"label"] = 
    }
    else{
      if (x[gene, "nongen_high"] == 1 ){
        x[gene,"color"] = "orange"}
      if (x[gene, "nongen_low"] == 1 ){
        x[gene,"color"] = "green"}
      if (x[gene, "gen_missense"] == 1 ){
        x[gene,"color"] = "purple"}
      if (x[gene, "gen_amplification"] == 1 ){
        x[gene,"color"] = "Blue"}
      if (x[gene, "gen_truncating"] == 1 ){
        x[gene,"color"] = "red"}
    }
  }
  ############################################################################
  ########################### label modifier ###########################
  ############################################################################
  
  gene = 1
  for (gene in genes){
    if (x[gene, "color"] == "black" ){
      postfix = ""
      if (x[gene, "nongen_high"] == 1){
        paste0(postfix, "#")
        # x[gene,"label"] = paste0(row.names(x[gene]," #"))
      }
      if (x[gene, "nongen_low"] == 1){
        postfix = paste0(postfix, "*")
        # x[gene,"label"] = paste0(row.names(x[gene]," *"))
      }
      if (x[gene, "gen_missense"] == 1){
        postfix = paste0(postfix, "@")
        # x[gene,"label"] = paste0(row.names(x[gene]," @"))
      }
      if (x[gene, "gen_amplification"] == 1){
        postfix = paste0(postfix, "$")
        # x[gene,"label"] = paste0(row.names(x[gene]," $"))
      }
      if (x[gene, "gen_truncating"] == 1){
        postfix = paste0(postfix, "%")
        # x[gene,"label"] = paste0(row.names(x[gene]," $"))
      }
      x[gene,"label"] = paste0(x[gene,"label"], postfix)
    }
  }
  x
}

plotInChromosomeContext <- function(candidGenes, geneAnnot, fname) {
  plot_title = fname
  
  gene_names = candidGenes
  
  genes_hg38_ = geneAnnot[which(geneAnnot$gene_name %in% rownames(gene_names)), c("seqnames", "start", "end", "gene_name", "width")]
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
  
  gene_names <- gene_names[order(match(row.names(gene_names),genes$gene_name)),]
  
  ############################################# even gene_names #############################################
  ids = seq(2, length(genes), by=2)
  tmp_genes = sapply(genes$gene_name[ids], toString)
  row_names_gene_names <- row.names(gene_names)
  colors2 <- vector()
  for (i in 1:length(row_names_gene_names)){
    for (j in 1:length(tmp_genes)){
      if (row_names_gene_names[i] == tmp_genes[j]){
        colors2 = c(colors2,gene_names$color[i])
      }}}
  
  for (i in 1:length(colors2)){
    if (colors2[i]=="black"){
      tmp_genes[i] = gene_names$label[i*2]
    }
  }
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 1, label.color=colors2, font=3)
  
  ############################################# odd gene_names #############################################
  ids = seq(1, length(genes), by=2)
  tmp_genes = sapply(genes$gene_name[ids], toString)
  row_names_gene_names <- row.names(gene_names)
  colors2 <- vector()
  for (i in 1:length(row_names_gene_names)){
    for (j in 1:length(tmp_genes)){
      if (row_names_gene_names[i] == tmp_genes[j]){
        colors2 = c(colors2,gene_names$color[i])
      }}}
  
  for (i in 1:length(colors2)){
    if (colors2[i]=="black"){
      tmp_genes[i] = gene_names$label[i*2-1]
    }
  }
  kpPlotMarkers(kp, data=genes[ids, ], labels=tmp_genes, r1 = 0.6, cex=0.65, text.orientation = "horizontal", data.panel = 2, label.color=colors2, font=3)
  legend("bottomright", legend=c("Non-genetic tumor suppersor", "Non-genetic oncogene", "Mutation oncogene","Amplification oncogene", "Mutation tumor suppersor", "Multiple potential drivers"), 
         col=c("orange", "green", "purple", "blue", "red", "Black"), cex=1, lty=1, lwd=2)
}

candidateGenes  <- function(cancer_type, adj_alpha=0.2 ) {
  
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-mutation-oncogene.csv"), stringsAsFactors = FALSE)
  
  nrOfCelllines = dat[1, ]$X.wt + dat[1,]$X.mut
  minNrOfCelllines=ceiling(nrOfCelllines/100)
  
  gene_set1 = gene_set2 = ""
  # non-genetic high : potential tumor suppressor
  fname = paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                 "-non-genetic-tumor_suppressor.csv")
  if(file.exists(fname)) {
    dat = read.csv(fname, stringsAsFactors = FALSE, row.names = 1)
    n = sum(dat$X.wt > minNrOfCelllines)
    # indexes = which(dat$p_adj < adj_alpha & dat$tcga_pvalue < 0.05 & dat$tcga_tumor_vs_normal == "down")
    indexes = which(dat$pvalue < 1/n & dat$tcga_pvalue < 0.05 & dat$tcga_tumor_vs_normal == "down")
    gene_set1 = row.names(dat)[indexes]
  }
  
  # non-genetic low
  fname = paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                 "-non-genetic-oncogene.csv")
  if(file.exists(fname)) {
    dat = read.csv( fname, stringsAsFactors = FALSE, row.names = 1)
    n = sum(dat$X.wt > minNrOfCelllines)
    # indexes = which(dat$p_adj < adj_alpha & dat$tcga_pvalue < 0.05 & dat$tcga_tumor_vs_normal == "up")
    indexes = which(dat$pvalue < 1/n & dat$tcga_pvalue < 0.05 & dat$tcga_tumor_vs_normal == "up")
    gene_set2 = row.names(dat)[indexes]
  }
  
  # missense
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-mutation-oncogene.csv"), stringsAsFactors = FALSE, row.names = 1)
  # indexes = which(dat$p_adj < adj_alpha)
  indexes = which(dat$pvalue < 1/n)
  gene_set3 = row.names(dat)[indexes]

  # amplification
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-amplification-oncogene.csv"), stringsAsFactors = FALSE, row.names = 1)
  # indexes = which(dat$p_adj < adj_alpha)
  indexes = which(dat$pvalue < 1/n)
  gene_set4 = row.names(dat)[indexes]
  
  # truncating
  dat = read.csv(paste0("../apsic_shiny/apsic_pvalues/", cancer_type, "/", cancer_type, 
                        "-mutation-tumor_suppressor.csv"), stringsAsFactors = FALSE, row.names = 1)
  # indexes = which(dat$p_adj < adj_alpha)
  indexes = which(dat$pvalue < 1/n)
  gene_set5 = row.names(dat)[indexes]

  list(nongen_high=gene_set1, nongen_low=gene_set2, gen_missense=gene_set3, 
       gen_amplification=gene_set4, gen_truncating=gene_set5)
}

#################################################################################################################################################

geneAnnot = read.table("Homo_sapiens.GRCh38.p10_genes.csv", sep='\t', stringsAsFactors = FALSE)

output_folder =  "figures/fig3/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
file_names = list.files("../apsic_shiny/apsic_pvalues/")
file_names <- file_names[-which(file_names == "Pan_cancer")]

nr_elem_candidate_genes <- function(candidGenes) {
  return(length(candidGenes$nongen_high) + length(candidGenes$nongen_low) + length(candidGenes$gen_missense) +
  length(candidGenes$gen_amplification) + length(candidGenes$gen_truncating))
}

for(fname in file_names){
  candidGenes = candidateGenes(fname)
  if(nr_elem_candidate_genes(candidGenes) > 0)  {
    gene = c(candidGenes$nongen_low,candidGenes$nongen_high,candidGenes$gen_missense,candidGenes$gen_amplification,candidGenes$gen_truncating)
    candidGenes <- annotateGene(candidGenes)
    pdf(paste0(output_folder, fname, ".pdf"), 10, 15)
    plotInChromosomeContext(candidGenes, geneAnnot, fname)
    dev.off()
  }
}

