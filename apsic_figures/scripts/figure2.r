rm(list=ls())
library(pvclust)
source("common_functions.r")

fig_folder = "figures/fig2/"
dir.create(fig_folder, recursive = T, showWarnings = FALSE)

##### we first load all p-value files
cancers = list.files("../apsic_shiny/apsic_pvalues/")
cancers = cancers[-which(cancers == "Pan_cancer")]

pvalues = prepareNongeneticPvalueData("non-genetic-tumor_suppressor", cancers)

################# clustering cancer subtypes according to the p-values ###############################
######################################################################################################
### first select 200 genes with highly variable p-values among subtypes
sorted_genes = sort(apply(pvalues, 1, var), decreasing = TRUE)
genes = names(sorted_genes)[1:500]
dat = pvalues[genes, ]

# select NA values to 1, if any. (In our case, there is no NA p-value so the following command has no impact)
dat[which(is.na(dat))] = 1.0

# clustering with 10000 bootstraps 
result = pvclust(dat, method.dist="cor", method.hclust="ward.D2", nboot=1000)

pdf(paste0(fig_folder, "types-cluster.pdf"), width=10, height=8)
plot(result)
pvrect(result, alpha=0.9)
dev.off()

# now we reorder subtypes according to the hierarchical clustering
clusters = hclust(get_dist(t(dat), method = "pearson"), method= "ward.D2")
plot(clusters)
cancers = clusters$labels[clusters$order]
save(cancers, file="ordered_cancers.RData")


################################## heatmap for 5 top genes per cancer ###############################
#####################################################################################################
rm(list=ls())
source("common_functions.r")
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

