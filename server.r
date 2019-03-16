
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
  output$p_wt_bladder_amp_low <- renderText ({
    dat = read.csv("apsic_pvalues/Bladder_Carcinoma/Bladder_Carcinoma-amplification-low.csv", row.names = 1)
    dat[input$gene, "pvalue_wt"]
  
  })
  
  output$p_wt_bladder_mis_low <- renderText ({
    dat = read.csv("apsic_pvalues/Bladder_Carcinoma/Bladder_Carcinoma-missense-low.csv", row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_wt_bladder_non_gen_hi <- renderText ({
    dat = read.csv("apsic_pvalues/Bladder_Carcinoma/Bladder_Carcinoma-non-genetic-high.csv", row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_wt_bladder_non_gen_low <- renderText ({
    dat = read.csv("apsic_pvalues/Bladder_Carcinoma/Bladder_Carcinoma-non-genetic-low.csv", row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  output$p_wt_bladder_trun_hi <- renderText ({
    dat = read.csv("apsic_pvalues/Bladder_Carcinoma/Bladder_Carcinoma-truncating-high.csv", row.names = 1)
    dat[input$gene, "pvalue_wt"]
    
  })
  
  
  
}



)