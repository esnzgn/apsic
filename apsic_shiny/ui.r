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
              
              tabPanel("Chromose Karyoplots",tags$label("tabPanel('Chromose Karyoplots',tags$img(src = paste0('/', '/figures/Karyotype', cancer_type))")),
              tabPanel("What we did in APSIC", tags$video(src="genomics_panel.mp4", width="800px", height="600px", type="video/mp4", controls = "controls"),
                                               tags$video(src="Tehran_khabar.mp4", width="800px", height="600px", type="video/mp4", controls = "controls")),
              tabPanel("APSIC Algorithm", tags$img(src = "apsic_algo.png",width="700px",height="800px" )),
              tabPanel("About us",textOutput("about_us"))
              )
              
            )
              )
                  )
        )

# shiny_www_address in Ehsan's system    D:\Program Files\R\R-3.5.3\library\shiny\www