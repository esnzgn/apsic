rm(list=ls())
library(tidyverse)
library(grid)

# load viability data
load("../apsic_shiny/cancerData.RData")
celllinesTable = table(cancerData$CellLine_annot$PATHOLOGIST_ANNOTATION)
tab = celllinesTable
tab = celllinesTable[celllinesTable>=5]

tab = sort(tab, decreasing = TRUE)

value = as.numeric(c(unname(tab)))
individual=names(tab)
id =1:length(value)

data = data.frame(id=id, individual = individual, value = value)

p = ggplot(data, aes(reorder(individual, -value), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", fill=alpha("darkorange3", 0.7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  labs(x = "", y = "Frequency") +
  theme(plot.margin=unit(c(0.5,0,0,1),"cm"))

p
 
fig_folder = "figures/fig_supp/"
dir.create(fig_folder, recursive = T, showWarnings = FALSE)

pdf(paste0(fig_folder, "barplot_cancers.pdf"), 10, 6)
print(p)
dev.off()

png(paste0(fig_folder, "barplot_cancers.png"), width = 10, height = 6, units = 'in', res = 300)
print(p)
dev.off()
# 
# 
