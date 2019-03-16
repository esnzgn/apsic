rm(list=ls())
source("../apsic_shiny/common.r")
source("../apsic_shiny/waterfall_plot_methods.r")

# load viability data
load("../apsic_shiny/cancerData.RData")

set.seed(1397)

fig_folder = "figures/fig1/"
dir.create(fig_folder, recursive = T, showWarnings = FALSE)

gene = "TP53"

############## Figure 1-a
pdf(paste0(fig_folder, "fig-1-a.pdf"), 8, 6) 
waterfallForGene(cancerData, gene = gene, title="", rank=FALSE, legenedPos="bottomleft", 
                 cols=NULL, type="all", sig_alpha = NA)
dev.off()

# png(paste0(fig_folder,"fig-1-a.png", res = 300))
# par(bg=NA)
# waterfallForGene(cancerData, gene = gene, title="", rank=FALSE, legenedPos="bottomleft", 
#                  cols=NULL, type="all", sig_alpha = NA)
# dev.off()


############## Figure 1-b
pdf(paste0(fig_folder,"fig-1-b.pdf"), 8, 6)
waterfallForGene(cancerData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="all", sig_alpha = NA)
dev.off()

############## Figure 1-c
pdf(paste0(fig_folder,"fig-1-c.pdf"), 8, 6)
gene =3000
waterfallForGene(cancerData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="all", sig_alpha = NA)
dev.off()

############## Figure 1-d  -- KRAS
pdf(paste0(fig_folder,"fig-1-d_wt_kras.pdf"), 10, 6)
gene = "KRAS"
waterfallForGene(cancerData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_wt", sig_alpha = NA)

dev.off()

pdf(paste0(fig_folder, "fig-1-d_missense_kras.pdf"), 8, 6)
waterfallForGene(cancerData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_missense", sig_alpha = NA)
dev.off()

############## Figure 1-d  -- TP53
pdf(paste0(fig_folder,"fig-1-d_wt_tp53.pdf"), 5, 4)
gene = "TP53"
waterfallForGene(cancerData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_wt", sig_alpha = NA)

dev.off()

pdf(paste0(fig_folder,"fig-1-d_missense_tp53.pdf"), 5, 4)
waterfallForGene(cancerData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_missense", sig_alpha = NA)
dev.off()


selectedData = selectCelllines(cancerData, "Breast:Carcinoma", tableS2File="../apsic_shiny/TableS2.csv")
gene = "NR3C2"
pdf(paste0(fig_folder,"fig-1-non-genetic-wt.pdf"), 5, 4)
waterfallForGene(selectedData, gene = gene, title="", rank=TRUE, legenedPos="bottomleft", 
                 cols=NULL, type="only_wt", sig_alpha = NA)
dev.off()



PARP4

