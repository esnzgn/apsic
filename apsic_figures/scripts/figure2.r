

################################## heatmap for 5 top genes per cancer ###############################
#####################################################################################################
rm(list=ls())
source("scripts/common_functions.r")
load("ordered_cancers.RData")

fig_folder = "figures/fig2/"
feature_type = "non-genetic-tumor_suppressor"
pvalues = prepareNongeneticPvalueData(feature_type, cancers)
topGenes = selectImportantNongeneticGenes(feature_type, cancers, 5, FALSE)
pvalues = pvalues[topGenes, ]

pdf(paste0(fig_folder, feature_type, ".pdf"), 8, 14)
plotHeatmap(-log10(pvalues) )
dev.off()

png(paste0(fig_folder, feature_type, ".png"), width = 8, height = 12, units = 'in', res = 300)
plotHeatmap(-log10(pvalues) )
dev.off()


feature_type = "non-genetic-oncogene"
pvalues = prepareNongeneticPvalueData(feature_type, cancers)
topGenes = selectImportantNongeneticGenes(feature_type, cancers, 5, FALSE)
pvalues = pvalues[topGenes, ]

pdf(paste0(fig_folder, feature_type, ".pdf"), 8, 14)
plotHeatmap(-log10(pvalues) )
dev.off()

png(paste0(fig_folder, feature_type, ".png"), width = 8, height = 12, units = 'in', res = 300)
plotHeatmap(-log10(pvalues) )
dev.off()

