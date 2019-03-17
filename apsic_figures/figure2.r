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
# dev.off()

# pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE)
# 
# pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE,
#          color = colorRampPalette(c("white","yellow2"))(10), 
#          border_color = TRUE)
# 
# pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE,
#          color = colorRampPalette(c("white","yellow2"))(10), 
#          border_color = "yellow2")
# 
# pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE,
#         color = colorRampPalette(c("white","yellow2"))(10), 
#         border_color = "white")

pheatmap(dat, cluster_rows = FALSE,cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(10), 
         border_color = FALSE, display_numbers = matrix(ifelse(dat > 0.95, "5%", ""), nrow(dat)),
         legend = TRUE, angle_col = "45")