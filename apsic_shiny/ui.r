shinyUI(fluidPage(
  titlePanel("APSiC Data Portal"),
  sidebarLayout(
    sidebarPanel(("Enter gene/cancer information"),
                 selectInput("gene",   "please select gene symbol", choices = rownames(cancerData$viabilities)),
                 selectInput("cancer", "Please select cancer type", choices = c("Pan cancer", cancer_types)),
                 radioButtons("filter", label = "Select the filter", 
                              choiceNames  = list("Mutation", "Copy number"), 
                              choiceValues = list("mutation", "CNA"))),
    
    mainPanel(("Gene/Cancer chart Information"),
              
              h4("P-Values"),
              h5("p-value wild type cancer amplification-low"),
              textOutput("p_wt_cancer_amp_low"),
              
              h5("p-value wild type cancer missense-low"),
              textOutput("p_wt_cancer_mis_low"),
              
              h5("p-value wild type cancer non-genetic-high"),
              textOutput("p_wt_cancer_non_gen_hi"),
              
              h5("p-value wild type cancer non-genetic-low"),
              textOutput("p_wt_cancer_non_gen_low"),
              
              h5("p-value wild type cancer truncating-high"),
              textOutput("p_wt_cancer_trun_hi"),
              
              h5("p-value gene expression TCGA"),
              textOutput("p_TCGA"),
              
              plotOutput("barChart"),
              plotOutput("wfplot_Mut_CNV"),
              plotOutput("wfplot_only_mut"),
              plotOutput("wfplot_only_wt")
              
            
              
              
    )
    
  )
  
))