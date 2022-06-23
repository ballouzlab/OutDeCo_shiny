server <- function(input, output) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })

  # renders table to review uploaded data in Assess DE
#  output$contents <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
#    inFile <- input$file1
    
#    if (is.null(inFile))
#      return(NULL)
    
#    read.csv(inFile$datapath)
#  })

  data <- reactive({
  inFile <- input$file2
  if(is.null(inFile)){return()} 
  read.table(file=inFile$datapath, sep=input$sep)
  })

  output$contents <- renderTable({
  if(is.null(data())){return ()}
  data()
  })

  output$tb <- renderUI({
    tableOutput("contents")
  })
















  # renders table to review uploaded data in Run DE
#  output$contents <- renderTable({
    # input$file2 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
#    inFile <- input$file2
#    ext <- tools::file_ext(inFile$datapath)
#    req(inFile)
#    validate(need(ext == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
#    read.csv(inFile$datapath)
#  })
}
