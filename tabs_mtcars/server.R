library(shiny)
shinyServer(function(input,output){
  output$mtcars <- renderTable({
    mtcars[,c("mpg",input$ngear)]
  })
  
  output$summ <- renderPrint({
    summary(mtcars[,c("mpg",input$ngear)])
  })
  
  output$plot <- renderPlot({
    with(mtcars,boxplot(mpg:gear))
  })
})