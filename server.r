rm(list=ls())
library(shiny)
source("waterfall_plot_methods.r")
source("common.r")
folder = ""
load(paste0(folder, "cancerData.RData"))



shinyServer(function(input, output){
  output$GS <- renderText(input$gene)
  output$CN <- renderText(input$cancer)
  output$filter <- renderText(input$filter)
  
 
 
  
  
  output$rank <- renderText(input$rank)
  output$wf0plot <- renderPlot (
    
    waterfallForGene(cancerData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft", 
                     cols=NULL, type = input$filter, sig_alpha = NA))
  
  output$wfplot <- renderPlot ( {
    # load tcga data
    cancer_type = input$cancer
    cancerToTCGA_data = read.csv("CelllinesToTumorTypes.csv", stringsAsFactors = FALSE)
    rownames(cancerToTCGA_data) = cancerToTCGA_data$Celllines
    tmp = cancerToTCGA_data[cancer_type, "mappedTCGAFile"]
    if(is.na(tmp) == FALSE) {
      load(paste0("tcga/", tmp, ".RData"))
      boxplot_gene(tcga_data, input$gene)
    }
  }
    )
  
  
}



)
