source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)

server <- function(input, output, session) {
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })
  
  observe({
      # DEFile from fileInput() function
      ServerDEFile <- req(input$DEFile)
      
      # extensions tool for format validation
      extDEFile <- tools::file_ext(ServerDEFile$datapath)
      if (is.null(input$DEFile)) {
        return ()
      } else{
        if (extDEFile == "txt") {
          label = paste("Delimiters for", extDEFile, "file")
          choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
        } else if (extDEFile == "tsv") {
          label = paste("Delimiter: Tab")
          choice <- (Tab="\t")
        } else {
          label = paste("Delimiter: Comma")
          choice <- (Comma=",")
        }
        updateRadioButtons(session, "sepButton", label = label, choices = choice)
      }
    })
  
  # reactive converts the upload file into a reactive expression known as data
  DEData <- reactive({

    # DEFile from fileInput() function
    ServerDEFile <- input$DEFile

    # extensions tool for format validation
    extDEFile <- tools::file_ext(ServerDEFile$datapath)

    # file format checking
    req(ServerDEFile)
     validate(need(extDEFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))

    # convert data into file format
    if (is.null(extDEFile)) {
      return ()
    }

    read.table(file=ServerDEFile$datapath, sep=input$sepButton, header=TRUE, nrows=5)
  })

  # creates reactive table called DEFileContent
  output$DEFileContent <- renderTable({
    if (is.null(DEData())) {
      return ()
    }
    DEData()
  })

  # handles rendering of reactive object on tb on ui
  output$UIDEContent <- renderUI({
    tableOutput("DEFileContent")
  })

  # load gene_list 
  gene_list <- reactive(
    if (input$gene_list_selection == "Generate Gene List") {
      sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo )
    } else {
      read.delim(file = input$DEFile$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
    }
  )

  # generate subnetwork
  observeEvent(
    input$generate_subnet, 
    {

      output$subnetwork <- renderTable({
        sub_nets <- subset_network_hdf5_gene_list(gene_list(), tolower(input$network_type), dir="../networks/")
      })

      s = input$subnetwork_search
      txt = if (is.null(s) || s == '') 'Filtered data' else {
        sprintf('Data matching "%s"', s)
      }
    },
  )
  sn <- reactiveValues(sub_nets = NULL)
  # cluster genes
  observeEvent(
    {input$run},
    {
    
    # # 1. Generate/load gene list eg. X chromosome genes 
    # gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr=="chrX"], 100 )

    # # 2. Extract sub-network and other network properties
    # network_type <- 'generic'
    # sub_nets <- subset_network_hdf5_gene_list(gene_list, network_type, dir="../networks/")

    # # 3. Extract results
    # sub_net <- sub_nets$sub_net
    # node_degrees <- sub_nets$node_degrees 
    # medK <- as.numeric(sub_nets$median)

    # # 4. Run clustering  
    # clust_net <- list() 
    # clust_net[["genes"]] <- cluster_coexp( sub_net$genes, medK = medK, flag_plot = FALSE )

    sn$sub_nets <- subset_network_hdf5_gene_list(gene_list(), tolower(input$network_type), dir="../networks/")

    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees
    medK <- as.numeric(sn$sub_nets$median)

    clust_net <- list() 
    clust_net[["genes"]] <- cluster_coexp( sub_net$genes, medK = medK, flag_plot = FALSE )

    # add heading here 
    output$CNtext = renderText({
      input$run
      req(input$run) #to prevent print at first lauch
      isolate(print("Network of Clustered Genes"))
    }) 
    
    # network output
    output$network <- renderPlot(
      {plot_network( sub_net$gene, clust_net$gene, medK )},
      # width = 500,
      # height = 500
    )

    output$CHtext = renderText({
      input$run
      req(input$run) #to prevent print at first lauch
      isolate(print("Heatmap of Clustered Genes"))
    })


    # heatmap output
    output$heatmap <- renderPlot(
      {plot_coexpression_heatmap( sub_net$gene, clust_net$gene, flag_plot_bin = FALSE)},
      width = 500,
      height = 500
    )

    output$CHBtext = renderText({
      input$run
      req(input$run) #to prevent print at first lauch
      isolate(print("Binarized Heatmap of Clustered Genes"))
    })


    # heatmap output
    output$Bheatmap <- renderPlot(
      {plot_coexpression_heatmap( sub_net$gene, clust_net$gene )},
      width = 500,
      height = 500
    )

  })
  
  # GENE CONNECTIVITY
   observeEvent(
    {input$runGC},
    {

    # Run clustering if not done previously
    if (is.null(sn$sub_nets)) {
      sn$sub_nets <- subset_network_hdf5_gene_list(gene_list(), tolower(input$network_type), dir="../networks/")
        
    }

    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees  

    # Text output for density plot
    output$GCdensityGtext = renderText({
      input$runGC
      req(input$runGC) #to prevent print at first lauch
      isolate(print("Density plot of Gene Connectivity"))
    }) 

    # density output
    output$GCdensityG <- renderPlot(
      {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                  node_degrees$genes[,2]/node_degrees$n_genes, 
                  xlab="Global node degree", 
                  ylab="Local node degree", flag= "density")  },
       width = 500,
       height = 500
    )

    # Text output for density plot
    output$GChistogramGtext = renderText({
      input$runGC
      req(input$runGC) #to prevent print at first lauch
      isolate(print("Histogram of Gene Connectivity"))
    }) 

    # histogram output
    output$GChistogramG <- renderPlot(
      {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                  node_degrees$genes[,2]/node_degrees$n_genes, 
                  xlab="Global node degree", 
                  ylab="Local node degree", flag= "hist")  },
       width = 500,
       height = 500
    )

    # m <- match(clust_net$genes$clusters$genes , rownames(sub_net$genes))

    # output$GChistogramSubsetG <- renderPlot(
    # { plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
    #                node_degrees$genes[m,2]/node_degrees$n_genes, 
    #                xlab="Global node degree", 
    #                ylab="Local node degree", 
    #                clusters = clust_net$genes$clusters )   },
    #     width = 500,
    #     height = 500
    # )



  })
}


