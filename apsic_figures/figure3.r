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
  legend("bottomright", legend=c("nongen_high #", "nongen_low *", "gen_missense @","gen_amplification $", "gen_truncating %", "multiple modifications: #,*,@,$,%"), 
         col=c("orange", "green", "purple", "blue", "red", "Black"), cex=1, lty=1, lwd=2)
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

#################################################################################################################################################

geneAnnot = read.table("Homo_sapiens.GRCh38.p10_genes.csv", sep='\t', stringsAsFactors = FALSE)

output_folder =  "figures/fig3/"
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
file_names = list.files("../apsic_shiny/apsic_pvalues/")
file_names <- file_names[-21]

for(fname in file_names){
  candidGenes = candidateGenes(fname)
  gene = c(candidGenes$nongen_low,candidGenes$nongen_high,candidGenes$gen_missense,candidGenes$gen_amplification,candidGenes$gen_truncating)
  candidGenes <- annotateGene(candidGenes)
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
