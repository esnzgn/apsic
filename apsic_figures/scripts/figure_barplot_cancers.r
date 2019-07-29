rm(list=ls())
library(tidyverse)


# load viability data
load("../apsic_shiny/cancerData.RData")
celllinesTable = table(cancerData$CellLine_annot$PATHOLOGIST_ANNOTATION)
tab = celllinesTable[celllinesTable>=5]

# load cancers
load("ordered_cancers.RData")

#### for debgging purpose
namesTable = cbind(sort(cancers), sort(names(tab)),  str_replace(sort(cancers), "_", ":") )
namesTable[which(namesTable[, 2] != namesTable[, 3])]
####

# changing "_" to ":"
cancers = str_replace(cancers, "_", ":")
# easier to handle these special case
cancers[which(cancers == "Upper:Aerodigestive_Tract_Carcinoma")] = "Upper_Aerodigestive_Tract:Carcinoma"
cancers[which(cancers == "Soft:Tissue_Sarcoma_Rhabdoid")] = "Soft_Tissue:Sarcoma_Rhabdoid"


tab = tab[cancers]

value = c(unname(tab))
individual=names(tab)
id =1:length(value)

data = data.frame(id=id, individual = individual, value = value)

barplot(data)

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data=data

# calculate the ANGLE of the labels
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
label_data$label<-ifelse(angle < -90, paste0(individual, "  ", value), paste0( value,  "  ", individual)) 
# ----- ------------------------------------------- ---- #


# Start the plot
p = ggplot(data, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-100,120) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=value+10, label=label_data$label, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p

fig_folder = "figures/fig_supp/"
dir.create(fig_folder, recursive = T, showWarnings = FALSE)


pdf(paste0(fig_folder, "circ_barplot_cancers.pdf"), 8, 8)
print(p)
dev.off()

png(paste0(fig_folder, "circ_barplot_cancers.png"), width = 8, height = 8, units = 'in', res = 300)
print(p)
dev.off()
