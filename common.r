
IH_CDF <- function(x, n) {
  if(n <= 20){
    X <-  floor(x)
    k <- seq(from = 0, to = X)
    # compute the cdf or p-value for n <= 100
    s <-  (-1)^k * choose(n, k)*( (x-k)^n)
    return(sum(s)/factorial(n))
  }else{
    # approximation for large n
    return(pnorm(x, n/2,sqrt(n/12), lower.tail = TRUE))
  }
}

getCancerTypes <- function() {
  base_folder = "Data/"
  meta_data <- read.csv(paste0(base_folder, "MutationFiles/File_metadata.csv"), header = TRUE, stringsAsFactors = FALSE)
  meta_data$Primary_site
}




# selectCelllines <- function(panCancerData, text, filter=NA) {
selectCelllines <- function(panCancerData, patho_annot) {
  CellLine_annot <-  read.csv("Data/ProjectDRIVE/TableS2.csv",stringsAsFactors = F)
  CellLine_annot = CellLine_annot[CellLine_annot$PATHOLOGIST_ANNOTATION == patho_annot, , drop=FALSE]
  
  if(nrow(CellLine_annot) == 0 ) {
    return(NA)
  }
  
  # if(is.na(filter)) {
  #   CellLine_annot = CellLine_annot[CellLine_annot$PRIMARY_SITE == text , ]
  #   
  # } else {
  #   CellLine_annot = CellLine_annot[CellLine_annot$PRIMARY_SITE == text & CellLine_annot$PATHOLOGIST_ANNOTATION == filter, ]
  # }
  
  
  
  indexes = which(colnames(panCancerData$viabilities) %in% paste0(CellLine_annot$CELLLINE, "_", CellLine_annot$PRIMARY_SITE))
  
  panCancerData$viabilities = panCancerData$viabilities[, indexes, drop=FALSE]
  panCancerData$mutations_all = panCancerData$mutations_all[, indexes, drop=FALSE]
  panCancerData$silentMutations = panCancerData$silentMutations[, indexes, drop=FALSE]
  panCancerData$missenseMutations = panCancerData$missenseMutations[, indexes, drop=FALSE]
  panCancerData$truncatingMutations = panCancerData$truncatingMutations[, indexes, drop=FALSE]
  
  panCancerData$copyNumbers = panCancerData$copyNumbers[, indexes, drop=FALSE]
  # panCancerData$exprData = panCancerData$exprData[, indexes]
  panCancerData
}


##### plotting tcga data
##### plotting tcga data
boxplot_gene <- function(tcga_data, gene) {
  gene_index = which(rownames(tcga_data$normal_counts)== gene) 
  
  tumor_cpm = cpm(tcga_data$tumor_counts)
  normal_cpm = cpm(tcga_data$normal_counts)
  
  N = normal_cpm[gene_index, ] + 1
  T = tumor_cpm[gene_index, ] + 1
  
  dat = data.frame(group = c(rep( "Normal", length(N)), rep( "Tumor", length(T)) ) , CPM = c(N, T) ) 
  
  # geom_boxplot proposes several arguments to custom appearance
  gPlot = ggplot(dat, aes(x=group, y=CPM)) + 
    geom_boxplot(
      
      # custom boxes
      color="blue",
      fill="blue",
      alpha=0.2,
      
      # Notch?
      # notch=TRUE,
      # notchwidth = 0.8,
      
      # custom outliers
      outlier.colour="red",
      outlier.fill="red",
      outlier.size=3
      
    ) + scale_y_continuous(trans='log10') + theme(axis.title.x=element_blank())
  
  print(gPlot)
}


mapPrimarySiteToTCGAproject <- function(primary_site) {
  data = read.csv("Data/TCGA/TCGAProjectsDesc.csv", stringsAsFactors = F)
  rownames(data) = data$primary_site

  data[primary_site, ]$tcga_project
}
