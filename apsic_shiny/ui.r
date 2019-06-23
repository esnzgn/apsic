shinyUI(fluidPage(
  titlePanel("APSiC Data Portal"),
  sidebarLayout(
    sidebarPanel(fluidRow(column(12,
                                 selectInput("gene",   "please select gene symbol", choices = rownames(cancerData$viabilities),selected = "TP53"),
                                 selectInput("cancer", "Please select cancer type", choices = c("Pan cancer", cancer_types),selected = "Breast:Carcinoma"),
                                 radioButtons("filter", label = "Select the filter", 
                                              choiceNames  = list("Mutation", "Copy number"), 
                                              choiceValues = list("mutation", "CNA")))),
                 fluidRow(imageOutput("imageBody",width = "80%", height = "100%"),class="center"), width=3),
    
    
    mainPanel(tabsetPanel(type = "tab",
                          tabPanel("Water-fall plots",
                                   fluidRow(
                                     column(6, plotOutput("wfplot_only_wt")),
                                     column(6, plotOutput("wfplot_Mut_CNV"))), 
                                   fluidRow(
                                     column(4, plotOutput("wfplot_only_mut")),
                                     column(4, plotOutput("barChart")), 
                                     column( dataTableOutput("ptable"), width= 4,height= 4))),
                          
                          # tabPanel("Chromose Karyoplots",tags$label("tabPanel('Chromose Karyoplots',tags$img(src = paste0('/', '/figures/Karyotype', cancer_type))")),
                          # navbarMenu("APSIC",
                          #            tabPanel("What we did in APSIC", tags$video(src="sample_video1.mp4", width="800px", height="600px", type="video/mp4", controls = "controls"),
                          #                     tags$video(src="sample_video.mp4", width="800px", height="600px", type="video/mp4", controls = "controls")),
                                     tabPanel("About APSIC",class="text-center", tags$img(src = "apsic_algo.png",width="100%",height="100%"))
                                     # )
                          # navbarMenu("APSICian",
                          #            tabPanel("About us",textOutput("about_us")),
                          #            tabPanel(class="text-center","Contact us",fluidRow(tags$label("/r /n Email, Website, Citation, Social Media ...")),fluidRow(tags$img(src= "SM.png")))
                          #            )
              )
              
    )
    
  )
)
)