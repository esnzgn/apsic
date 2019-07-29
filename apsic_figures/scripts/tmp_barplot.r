rm(list=ls())
library(ggplot2)
library(reshape2)
library(viridis)
source("../apsic_shiny/common.r")
source("../apsic_shiny/waterfall_plot_methods.r")
source("scripts/plot_functions_for_profiles.r")
source("scripts/common_functions.r")

# load viability data
load("../apsic_shiny/cancerData.RData")
load("ordered_cancers.RData")

feature_type = "non-genetic-tumor_suppressor"
pvalues = prepareNongeneticPvalueData(feature_type, cancers)

head(pvalues)
pdf("a.pdf", width=15, height=15)

par(mar=c(10,4,4,2))
AA = apply(pvalues, 2, function(x) { sum(x < 0.05, na.rm = T) })
barplot(AA , border="white", names.arg = colnames(pvalues),   beside=T, las=2, cex.names=0.7)
dev.off()
#Create data
set.seed(112)
data=matrix(sample(1:30,15) , nrow=3)
colnames(data)=c("A","B","C","D","E")
rownames(data)=c("var1","var2","var3")

# Grouped barplot
barplot(data, col=colors()[c(23,89,12)] , border="white", font.axis=2, beside=T, legend=rownames(data), xlab="group", font.lab=2)
