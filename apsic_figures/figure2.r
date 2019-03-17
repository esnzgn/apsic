rm(list=ls())
library(pheatmap)

set.seed(10)

prepareData <-function(r, c) {
  data = matrix(runif(n=r*c), r, c)
  
  rownames(data) =paste0(rep(c("KRAS", "TP53", "MDM2", "MDM4", "GATA3"), r/5), 1:c)
  colnames(data)=paste0("TCGA", 1:c)
  data
}




# plotHeatmap <- function(data) {
#   # heat
#   
# }


dat = prepareData(20, 20)

# dat = matrix(rnorm(50*17),50,17)
# heatmap(dat)
dev.off()

# pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE)
pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE,
        color = colorRampPalette(c("black","firebrick3","white"))(10), 
        display_numbers = matrix(ifelse(dat > 0.8, "*", ""), nrow(dat)), border_color = FALSE)
