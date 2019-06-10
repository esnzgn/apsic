library(shiny)
shinyUI(fluidPage(                                  #a kind of page
  headerPanel(title = "Shiny Tabset Example"),
  sidebarLayout(                                    #a kind of page layout
    sidebarPanel(
      selectInput("ngear", "select the gear number",c("cylinders" = "cyl", "Transmission" = "am", "Gears" = "gear"))
      
    ),
    mainPanel(
      tabsetPanel(type = "tab",
                  tabPanel("Help",tags$img(src = "1.png")),
                  tabPanel("Data", tableOutput("mtcars")),
                  tabPanel("Summary",verbatimTextOutput("summ")),
                  tabPanel("plot",plotOutput("plot"))
                  )
    )
  )
))
