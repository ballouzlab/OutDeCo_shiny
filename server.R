server <- function(input, output) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })

  # reactive converts the upload file into a dataframe known as data
  data <- reactive({

  # DEFile from fileInput() function
  ServerDEFile <- input$DEFile

  # extensions tool for format validation
  ext <- tools::file_ext(ServerDEFile$datapath)

  # file format checking
  req(ServerDEFile)
  validate(need(ext == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

  # convert data into file format
  if(is.null(ServerDEFile)){return()} 
  read.table(file=ServerDEFile$datapath, sep=input$sep)
  })

  # creates reactive table called DEFileContent
  output$DEFileContent <- renderTable({
  if(is.null(data())){return ()}
  data()
  })

  # handles rendering of reactive object on tb on ui
  output$UIDEContent <- renderUI({
    tableOutput("DEFileContent")
  })
}
