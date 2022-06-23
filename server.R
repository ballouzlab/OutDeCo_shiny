server <- function(input, output) {
  theme = bs_theme(version = 4, bootswatch = "minty")
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })
  output$selected_var <- renderText({ 
    "You have selected this"
  })
  output$table <- DT::renderDT({
    iris
  })
  output$text <-  renderText({ output$txt })
}