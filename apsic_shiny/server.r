#install.packages("readr")
#install.packages("magick")
#install.packages("jpeg")
#install.packages("pdftools")
#install.packages("readbitmap")
# install.packages("ggplot2")
# install.packages("png")
# install.packages("DT")

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
      return(round(-log10(dat[gene, ptype]),2))
    }
    return(NULL)
  }
  
  output$ptable <- renderDataTable({ 
    
    
    
    data = data.frame(matrix(ncol=2, nrow=3))
    colnames(data) = c("Type", "-log10(p-value)")
    
    data[1, ] = c("amplification oncogene", getPValue(input$gene, input$cancer, "amplification-low",  "pvalue_mut") )
    data[2, ] = c("mutation oncogene", getPValue(input$gene, input$cancer, "missense-low",  "pvalue_mut")) 
    data[3, ] = c("mutation tumor suppressor", getPValue(input$gene, input$cancer, "truncating-high",  "pvalue_mut")) 
    
    p_value = getPValue(input$gene, input$cancer, "non-genetic-low",  "pvalue_wt")
    if(is.null(p_value) == FALSE) {
      data[4, ] = c("non-genetic tumor suppressor", p_value)    
      
    }
    p_value = getPValue(input$gene, input$cancer, "non-genetic-high",  "pvalue_wt")
    if(is.null(p_value) == FALSE) {
      data[5, ] = c("non-genetic oncogene", p_value)    
    }
    print(round(getPValue(input$gene, input$cancer, "non-genetic-low",  "pvalue_wt")),3)
    
    data}
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

