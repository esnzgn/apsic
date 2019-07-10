#install.packages("readr")
#install.packages("magick")
#install.packages("jpeg")
#install.packages("pdftools")
#install.packages("readbitmap")
# install.packages("ggplot2")
# install.packages("png")
# install.packages("DT")
# install.packages("shinydashboard")
# install.packages("shinyWidgets")

# library(jpeg)
# library(readbitmap)
# library(magick)
# # library(magrittr)
# library(imager)
# #library(pdftools)
library(readr)
# library("ggplot2")
library(shiny)
library(png)
library(DT)
# library(shinydashboard)
# library(shinyWidgets)


shinyServer(function(input, output,session){
  
  # output$addressBody <- renderText({
  #   if(input$cancer == "Breast:Carcinoma") {
  #     addressB = "/body/breast.png"
  #   } 
  #   if(input$cancer == "Pan cancer") {
  #     addressB = "/body/body.png"
  #   }
  #   if(input$cancer == "Leukemia:ALL") {
  #     addressB = "/body/blood.png"
  #   } else (addressB = "/body/body.png")
  #   addressB
  # })
  

  # waterfall_plot mutation or copy number
  
  output$wfplot_Mut_CNV <- renderPlot ({
    if(input$cancer == "Pan cancer") {
      selectedData = cancerData
      
    } else {
      selectedData = selectCelllines(cancerData, input$cancer)
    }
    
    if(input$filter == "mutation") {
      waterfallForGene(selectedData, gene = input$gene, title=paste0(input$gene,": all cell lines"), rank=TRUE, legenedPos="bottomleft",
                       cols=NULL, type = "all", sig_alpha = NA)
    } else if (input$filter == "CNA") {
      waterfallForGene_CNA(selectedData, gene = input$gene, title=paste0(input$gene,": all cell lines"), rank=TRUE, legenedPos="bottomleft",
                           cols=NULL, type = "all", sig_alpha = NA)
    }
  })
  
  # waterfall_plot only_mut
  output$wfplot_only_mut <- renderPlot ( {
    if(input$cancer == "Pan cancer") {
      selectedData = cancerData
      
    } else {
      selectedData = selectCelllines(cancerData, input$cancer)
    }
    
    if(input$filter == "mutation") {
      waterfallForGene(selectedData, gene = input$gene, title=paste0(input$gene,": altered cell lines"), rank=TRUE, legenedPos="bottomleft",
                       cols=NULL, type = "only_mut", sig_alpha = NA)
    } else if (input$filter == "CNA") {
      waterfallForGene_CNA(selectedData, gene = input$gene, title=paste0(input$gene,": altered cell lines"), rank=TRUE, legenedPos="bottomleft",
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
      waterfallForGene(selectedData, gene = input$gene, title=paste0(input$gene,": wild type cell lines"), rank=TRUE, legenedPos="bottomleft",
                       cols=NULL, type = "only_wt", sig_alpha = NA)
    } else if (input$filter == "CNA") {
      waterfallForGene_CNA(selectedData, gene = input$gene, title=paste0(input$gene,": wild type cell lines"), rank=TRUE, legenedPos="bottomleft",
                           cols=NULL, type = "only_wt", sig_alpha = NA)
    }
  }
  )
  
  
  output$barChart <- renderPlot ({
    # load tcga data barChart
    
    
    cancer_type = input$cancer
    cancerToTCGA_data = read.csv("CelllinesToTumorTypes.csv", stringsAsFactors = FALSE)
    rownames(cancerToTCGA_data) = cancerToTCGA_data$Celllines
    tmp = cancerToTCGA_data[cancer_type, "mappedTCGAFile"]
    if(is.na(tmp) == FALSE) {
      load(paste0("tcga/", tmp, ".RData"))
      boxplot_gene(tcga_data, input$gene)
    }
    ##########################################################################################################################################
    ##########################################################################################################################################
    ##########################################################################################################################################
    ##########################################################################################################################################
    if (input$cancer =="pan cancer"){
      return(NULL)
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
      return(round(-log10(dat[gene, ptype]),2))
    }
    return(NA)
  }
  
  getTCGAPValue <- function(gene, cancer) {
    tmp = str_replace(cancer, "[: ]", "_")
    # all p-value files have the same TCGA p-values 
    A <- paste0("apsic_pvalues/", tmp,"/", tmp, "-amplification-low.csv")
    if(file.exists(A)) {
      dat = read.csv(A, row.names = 1)
      pvalue = NA
      if("tumor_expressed_less" %in% colnames(dat)) {
        pvalue_less = dat[gene, "tumor_expressed_less"]
        pvalue_more = dat[gene, "tumor_expressed_more"]
        if(!is.na(pvalue_less)) 
          pvalue = round(-log10(min(pvalue_less, pvalue_more) ),2)
      }
      return(list(pvalue = pvalue))
    }
    return(NULL)
  }

  output$ptable <- DT::renderDataTable({ 
    print(input$cancer)
    data = data.frame(matrix(ncol=2, nrow=3))
    colnames(data) = c("Type", "-log10(p-value)")
    
    p_value = getPValue(input$gene, input$cancer, "amplification-low",  "pvalue_mut") 
    if(is.na(p_value))
      p_value = "-"
    
    data[1, ] = c("Amplification oncogene", p_value)
    
    p_value = getPValue(input$gene, input$cancer, "missense-low",  "pvalue_mut")
    if(is.na(p_value))
      p_value = "-"
    
    data[2, ] = c("Mutation oncogene", p_value) 
    
    p_value = getPValue(input$gene, input$cancer, "truncating-high",  "pvalue_mut")
    if(is.na(p_value))
      p_value = "-"
    
    data[3, ] = c("Mutation tumor suppressor", p_value) 
    
    if (input$cancer != "Pan cancer"){
      p_value = getPValue(input$gene, input$cancer, "non-genetic-low",  "pvalue_wt")
      if(is.na(p_value))
        p_value = "-"
      
      data[4, ] = c("Non-genetic oncogene", p_value)    
      
      p_value = getPValue(input$gene, input$cancer, "non-genetic-high",  "pvalue_wt")
      if(is.na(p_value))
        p_value = "-"
      
      data[5, ] = c("Non-genetic tumor suppressor", p_value)   
      
      res = getTCGAPValue(input$gene, input$cancer)
      p_value = res$pvalue 
      if(is.na(p_value))
        p_value = "-"
      
      data[6, ] = c("Comparison of TCGA expression levels for cancer and normal samples", p_value)   
    }
    
    return(data)
    }
    , options = list(dom = 't')
  )
  
  output$about_us <- renderText({
    # About us
    print(read_file("About_us.txt"))
  })
  
  output$imageBody <- renderImage({
    if (is.null(input$cancer))
      return(NULL)
    # print(input$cancer)
    
    if (input$cancer == "Breast:Carcinoma") {
      return(list(
        src = "body/breast.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Breast"
      ))}
    
     else if (input$cancer == "Leukemia:ALL") {
      return(list(
        src = "body/blood.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Leukemia"
      ))
    } else if (input$cancer == "CNS:Glioma_HighGrade") {
      return(list(
        src = "body/cns-pns.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "CNS"
      ))
    } else if (input$cancer == "Lung:NSCLC_Large_Cell") {
      return(list(
        src = "body/lung.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Lung"
      ))
    } else if (input$cancer == "Lung:NSCLC_Others") {
      return(list(
        src = "body/lung.png",width = "100%", height = "100%",
        contentType = "image/png",
        alt = "LUng_others"
      ))
    } else if (input$cancer == "Lung:SCLC") {
      return(list(
        src = "body/lung.png",width = "100%", height = "100%",
        contentType = "image/png",
        alt = "Lung_SC"
      ))
    } else if (input$cancer == "PNET:Neuroblastoma") {
      return(list(
        src = "body/cns-pns.png",width = "100%", height = "100%",
        contentType = "image/png",
        alt = "PNET:Neuroblastoma"
      ))
    } else if (input$cancer == "Soft_Tissue:Sarcoma_Rhabdoid") {
      return(list(
        src = "body/soft-tissue.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Soft_Tissue:Sarcoma_Rhabdoid"
      ))
    } else if (input$cancer ==  "Bladder:Carcinoma") {
      return(list(
        src = "body/Bladder.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Bladder"
      ))
    } else if (input$cancer ==  "Colorectal:Carcinoma") {
      return(list(
        src = "body/Colorectal.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Colorectal"
      ))
    } else if (input$cancer ==  "Lymphoma:NH_B_cell") {
      return(list(
        src = "body/Lymph.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Lymphoma:NH_B_cell"
      ))
    } else if (input$cancer ==  "Oesophagus:Carcinoma") {
      return(list(
        src = "body/Oesophagus.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Oesophagus:Carcinoma"
      ))
    } else if (input$cancer ==  "Upper_Aerodigestive_Tract:Carcinoma") {
      return(list(
        src = "body/aeroupperdigestive.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Upper_Aerodigestive_Tract:Carcinoma"
      ))
    } else if (input$cancer ==  "Kidney:Carcinoma") {
      return(list(
        src = "body/kidney.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Kidney:Carcinoma"
      ))
    } else if (input$cancer ==  "Leukemia:AML") {
      return(list(
        src = "body/blood.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Leukemia:AML"
      ))
    } else if (input$cancer == "CNS:Glioma") {
      return(list(
        src = "body/cns-pns.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "CNS:Glioma"
      ))
    } else if (input$cancer == "Liver:HCC") {
      return(list(
        src = "body/liver.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Liver:HCC"
      ))
    } else if (input$cancer =="Lung:NSCLC_Adeno") {
      return(list(
        src = "body/lung.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Liver:HCC"
      ))
    } else if (input$cancer == "Lung:NSCLC_Squamous") {
      return(list(
        src = "body/lung.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Lung:NSCLC_Squamous"
      ))
    } else if (input$cancer == "Lung:Mesothelioma") {
      return(list(
        src = "body/lung.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Lung:Mesothelioma"
      ))
    } else if (input$cancer == "Ovary:Carcinoma" ) {
      return(list(
        src = "body/ovary.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Ovary:Carcinoma"
      ))
    } else if (input$cancer == "Pancreas:Carcinoma" ) {
      return(list(
        src = "body/Pancreas.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Pancreas:Carcinoma"
      ))
    } else if (input$cancer == "Skin:Melanoma" ) {
      return(list(
        src = "body/skin.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Skin:Melanoma"
      ))
    } else if (input$cancer ==  "Gastric:Carcinoma") {
      return(list(
        src = "body/gastric.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Gastric:Carcinoma"
      ))
    } else if (input$cancer ==  "Thyroid:Carcinoma" ) {
      return(list(
        src = "body/Thyroid.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Thyroid:Carcinoma"
      ))
    } else if (input$cancer ==  "Endometrium:Carcinoma" ) {
      return(list(
        src = "body/Ovary.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Endometrium:Carcinoma"
      ))
    } else if (input$cancer ==  "Pan cancer" ) {
      return(list(
        src = "body/body.png",
        contentType = "image/png",width = "100%", height = "100%",
        alt = "Pan cancer"
      ))
    } 
    
  }
  , deleteFile = F)
  
  
})

