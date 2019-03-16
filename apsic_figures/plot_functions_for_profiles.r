plotDriverGenes <- function(MutMat) {
  
  ggplot(melt(MutMat), aes(Var1, Var2)) + 
    geom_tile(aes(fill = value),alpha=0.75)+
    ylab("Cell Lines") + 
    xlab("Genes") +
    theme_bw()+
    theme(panel.grid.major= element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_line(size = 0.5, color="#525252"),
          # axis.text.y=element_text(angle=0, hjust=1, size=11,colour="#525252"),
          # axis.text.x=element_text(angle=50, hjust=1, size=11,colour="#525252"),
          axis.text.y=element_blank(),
          axis.text.x= element_blank(),
          axis.title.x=element_text(angle=0, size=13, face="bold",vjust=1,colour="#525252"),
          axis.title.y=element_text(angle=90, size=13, face="bold",vjust=0,colour="#525252"),
          legend.position = "none")+
    scale_fill_gradient(low="white",high="red")
}


plotPerturbationProfile <- function(perturbData) {
  q <- ggplot(melt(perturbData), aes(Var1, Var2)) + 
    geom_tile(aes(fill = value),alpha=0.75)+
    ylab("Cell Lines") + 
    xlab("Perturbed genes") +
    theme_bw()+
    labs(fill = "viability") +
    theme(panel.grid.major= element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_line(size = 0.5,color="#525252"),
          # axis.text.y=element_text(angle=0, hjust=1,size=9,colour="#525252"),
          # axis.text.x=element_text(angle=50, hjust=1,size=9,colour="#525252"),
          axis.text.y=element_blank(),
          axis.text.x= element_blank(),
          axis.title.x=element_text(angle=0,size=12,face="bold",vjust=1,colour="#525252"),
          axis.title.y=element_text(angle=90, size=11,face="bold",vjust=0,colour="#525252"),
          legend.position = "bottom")+
    # scale_fill_viridis(option="plasma",begin = 0,end = 1)
    scale_fill_gradient2(low = "red", high="royalblue", mid="white", midpoint = 0)
  print(q)
}