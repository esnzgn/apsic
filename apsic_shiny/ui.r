shinyUI(fluidPage(
  titlePanel("APSiC Data Portal"),
  sidebarLayout(
    sidebarPanel(("Enter gene/cancer information"),
                 selectInput("gene",   "please select gene symbol", choices = rownames(cancerData$viabilities)),
                 selectInput("cancer", "Please select cancer type", choices = c("Pan cancer", cancer_types)),
                 radioButtons("filter", label = "Select the filter", 
                              choiceNames  = list("Mutation", "Copy number"), 
                              choiceValues = list("mutation", "CNA"))),
    
    mainPanel("Gene-Cancer Association",
              tabsetPanel(type = "tab",
              tabPanel("Water-fall plots",
                fluidRow(
                  column(6, plotOutput("wfplot_only_wt")),
                  column(6, plotOutput("wfplot_Mut_CNV"))), 
                fluidRow(
                  column(4, plotOutput("wfplot_only_mut")),
                  column(4, plotOutput("barChart")), 
                  column(4, tableOutput("ptable")))),
              tabPanel("Chromose Karyoplots",plotOutput("karyoplot")),
              tabPanel("APSIC Algorithm", tags$img(src = "rstudio.png")),
              tabPanel("About us",textOutput("about_us"))
              )
              
     
  )
  
)))
