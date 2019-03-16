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
              h5("p-value wild type Bladder Carcinoma amplification-low"),
              textOutput("p_wt_bladder_amp_low"),
              h5("p-value wild type Bladder Carcinoma missense-low"),
              textOutput("p_wt_bladder_mis_low"),
              h5("p-value wild type Bladder Carcinoma non-genetic-high"),
              textOutput("p_wt_bladder_non_gen_hi"),
              h5("p-value wild type Bladder Carcinoma non-genetic-low"),
              textOutput("p_wt_bladder_non_gen_low"),
              h5("p-value wild type Bladder Carcinoma truncating-high"),
              textOutput("p_wt_bladder_trun_hi"),
              
              
              
              plotOutput("barChart"),
              plotOutput("wfplot_Mut_CNV"),
              plotOutput("wfplot_only_mut"),
              plotOutput("wfplot_only_wt")
              
            
              
              
    )
    
  )
  
))