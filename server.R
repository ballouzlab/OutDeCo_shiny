server <- function(input, output, session) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })

  # reactive converts the upload file into a reactive expression known as data
  uploadData <- reactive({

  # DEFile from fileInput() function
  ServerDEFile <- input$DEFile

  # extensions tool for format validation
  extDEFile <- tools::file_ext(ServerDEFile$datapath)

  # file format checking
  req(ServerDEFile)
  validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

  # convert data into file format
  if(is.null(extDEFile)){return()} 

  if (extDEFile == "txt") {
    choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
    updateRadioButtons(session, "sepButton",
                      label = paste("Delimiter Selection for", extDEFile),
                      choices = choice
    )
  }
  else if (extDEFile == "csv") {
    choice <-c(",")
    updateRadioButtons(session, "sepButton",
                      label = paste("Delimiter: Comma"),
                      choices = choice
    )
  }
  else {
    choice <-c("\t")
    updateRadioButtons(session, "sepButton",
                      label = paste("Delimiter: Tab"),
                      choices = choice
    )
  }

  read.table(file=ServerDEFile$datapath, sep=input$sepButton)
  })

  # creates reactive table called DEFileContent
  output$DEFileContent <- renderTable({
  if(is.null(uploadData())){return ()}
  uploadData()
  })

  # handles rendering of reactive object on tb on ui
  output$UIDEContent <- renderUI({
    tableOutput("DEFileContent")
  })
}
