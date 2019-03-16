
shinyServer(function(input, output){
  output$GS <- renderText(input$gene)
  output$CN <- renderText(input$cancer)
  output$filter <- renderText(input$filter)
  output$rank <- renderText(input$rank)
  
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
  

  
  
  
}



)