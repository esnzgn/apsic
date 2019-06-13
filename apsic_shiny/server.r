#install.packages("readr")
#install.packages("magick")
#install.packages("jpeg")
#install.packages("pdftools")
#install.packages("readbitmap")
# install.packages("ggplot2")

# library(jpeg)
# library(readbitmap)
# library(magick)
# # library(magrittr)
# library(imager)
# #library(pdftools)
library(readr)
# library("ggplot2")
library(shiny)

shinyServer(function(input, output){
  
  
  output$GS <- renderText(input$gene)
  output$CN <- renderText(input$cancer)
  output$filter <- renderText(input$filter)
  output$rank <- renderText(input$rank)
  
  # output$GS <- renderText("TP53")
  # output$CN <- renderText("Breast:Carcinoma")
  
  # waterfall_plot mutation or copy number
  output$wfplot_Mut_CNV <- renderPlot ( {
    if(input$cancer == "Pan cancer") {
      selectedData = cancerData
      
    } else {
      selectedData = selectCelllines(cancerData, input$cancer)
    }
    
    if(input$filter == "mutation") {
      waterfallForGene(selectedData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft",
                       cols=NULL, type = "all", sig_alpha = NA)
    } else if (input$filter == "CNA") {
      waterfallForGene_CNA(selectedData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft",
                           cols=NULL, type = "all", sig_alpha = NA)
    }
  }
  )
  
  
  # waterfall_plot only_mut
  output$wfplot_only_mut <- renderPlot ( {
    if(input$cancer == "Pan cancer") {
      selectedData = cancerData
      
    } else {
      selectedData = selectCelllines(cancerData, input$cancer)
    }
    
    if(input$filter == "mutation") {
      waterfallForGene(selectedData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft",
                       cols=NULL, type = "only_mut", sig_alpha = NA)
    } else if (input$filter == "CNA") {
      waterfallForGene_CNA(selectedData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft",
                           cols=NULL, type = "only_cna", sig_alpha = NA)
    }
  }
  )
  
  
  
  # waterfall_plot only_wt
  output$wfplot_only_wt <- renderPlot ( {
    if(input$cancer == "Pan cancer") {
      selectedData = cancerData
      
    } else {
      selectedData = selectCelllines(cancerData, input$cancer)
    }
    
    if(input$filter == "mutation") {
      waterfallForGene(selectedData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft",
                       cols=NULL, type = "only_wt", sig_alpha = NA)
    } else if (input$filter == "CNA") {
      waterfallForGene_CNA(selectedData, gene = input$gene, title=input$gene, rank=TRUE, legenedPos="bottomleft",
                           cols=NULL, type = "only_wt", sig_alpha = NA)
    }
  }
  )
  
  
  output$barChart <- renderPlot ( {
    # load tcga data barChart
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
  
  #p_values
  output$p_wt_cancer_amp_low <- renderText ({
    tmp = str_replace(input$cancer, "[: ]", "_")
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-amplification-low.csv")
    dat = read.csv(A, row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_wt_cancer_mis_low <- renderText ({
    tmp = str_replace(input$cancer, "[: ]", "_")
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-missense-low.csv")
    dat = read.csv(A, row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_wt_cancer_non_gen_hi <- renderText ({
    tmp = str_replace(input$cancer, "[: ]", "_")
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-non-genetic-high.csv")
    dat = read.csv(A, row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_wt_cancer_non_gen_low <- renderText ({
    tmp = str_replace(input$cancer, "[: ]", "_")
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-non-genetic-low.csv")
    dat = read.csv(A, row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  
  output$p_wt_cancer_trun_hi <- renderText ({
    tmp = str_replace(input$cancer, "[: ]", "_")
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-truncating-high.csv")
    dat = read.csv(A, row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_TCGA <- renderText ({
    cancer_type = input$cancer
    cancerToTCGA_data = read.csv("CelllinesToTumorTypes.csv", stringsAsFactors = FALSE)
    rownames(cancerToTCGA_data) = cancerToTCGA_data$Celllines
    tmp = cancerToTCGA_data[cancer_type, "mappedTCGAFile"]
    if(is.na(tmp) == FALSE) {
      load(paste0("tcga/", tmp, ".RData"))
      return(tcga_data$pvalues[input$gene, 1])
    }
  })
  
  getPValue <- function(gene, cancer, alter_type, ptype ) {
    tmp = str_replace(cancer, "[: ]", "_")
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-", alter_type, ".csv")
    if(file.exists(A)) {
      dat = read.csv(A, row.names = 1)
      return(dat[gene, ptype])
    }
    return(NULL)
  }
  
  output$ptable <- renderTable({ 
    
    
    
    data = data.frame(matrix(ncol=2, nrow=3))
    colnames(data) = c("Type", "P Value")
    
    data[1, ] = c("amplification-low", getPValue(input$gene, input$cancer, "amplification-low",  "pvalue_mut") )
    data[2, ] = c("missense-low", getPValue(input$gene, input$cancer, "missense-low",  "pvalue_mut")) 
    data[3, ] = c("truncating-high", getPValue(input$gene, input$cancer, "truncating-high",  "pvalue_mut")) 
    
    p_value = getPValue(input$gene, input$cancer, "non-genetic-low",  "pvalue_wt")
    if(is.null(p_value) == FALSE) {
      data[4, ] = c("non-genetic-low", p_value)    
      
    }
    p_value = getPValue(input$gene, input$cancer, "non-genetic-high",  "pvalue_wt")
    if(is.null(p_value) == FALSE) {
      data[5, ] = c("non-genetic-high", p_value)    
    }
    print(round(getPValue(input$gene, input$cancer, "non-genetic-low",  "pvalue_wt")),3)
    
    data
  })
  
  output$about_us <- renderText({
    # About us
    print(read_file("About_us.txt"))
  })
})

