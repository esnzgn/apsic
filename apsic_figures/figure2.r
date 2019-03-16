rm(list=ls())

set.seed(10)

prepareData <-function() {
  data = matrix(runif(n=50*17), 50, 17)
  
  rownames(data) =paste0(rep(c("KRAS", "TP53", "MDM2", "MDM4", "GATA3"), 10), 1:15)
  colnames(data)=paste0("TCGA", 1:17)
  data
}

  
plotHeatmap function <- (data) {
  
}
