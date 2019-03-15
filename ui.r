library(shiny)

shinyUI(fluidPage(
  titlePanel("APSiC Data Portal"),
  sidebarLayout(
    sidebarPanel(("Enter gene/cancer information"),
                 selectInput("gene",   "please select gene symbol", choices = rownames(cancerData$viabilities)),
                 selectInput("cancer", "Please select cancer type", choices = c("Pan cancer", cancer_types)),
                 radioButtons("filter", label = "Select the filter", 
                              choiceNames  = list("Mutation", "Copy number"), 
                              choiceValues = list("all", "all"))),
    
    mainPanel(("Gene/Cancer chart Information"),
              
              h5("p-values amplification-low:0.01"),h5("p-values missense-low:0.01"),h5("p-values genetic-high:0.01"),h5("p-values genetic-low:0.01"), h5("p-values truncating-high:0.01"),
              plotOutput("wfplot"),
              plotOutput("wf0plot")
              
              
    )
    
  )
  
))