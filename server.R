source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)

server <- function(input, output, session) {
  # Removing elements that are not functional without subnetwork
  # hide(id = "runGC")
  # hide(id = "run")
  hide(id = "GC_dropdown")
  hide(id = "CG_dropdown")
  hide(id = "FO_dropdown")

  
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

  # Add the Run buttons 
  observeEvent(
    input$generate_subnet,
    {
      # cluster genes
      show(id = "CG_dropdown")
      show(id = "run")
      hide(id = "CG_error")

      # gene connectivity
      show(id = "GC_dropdown")
      show(id = "runGC")
      hide(id = "GC_error")

      # functional outliers
      show(id = "FO_dropdown")
      show(id = "runFO")
      hide(id = "FO_error")

    }
  )













  # sub_nets
  sn <- reactiveValues(
    sub_nets = NULL,
  )

  # generate sub_nets
  observeEvent(
    input$generate_subnet,
    {
      # determine gene_list
      if (input$gene_list_selection == "Generate Gene List") {
        validate(
          need(input$chooseChrome, 'Please enter a Chromosme between 1-22, X, Y'),
          need(input$chooseGeneNo != "" && input$chooseGeneNo > 0, 'Please enter a valid Gene Number'),
        )
        gene_list <- sample(EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo)
      } else {
        gene_list <- read.delim(file = input$DEFile$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
      }
      # generate sub_nets & store under sn
      sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(input$network_type), dir="../networks/")
    }
  )


  # clust_net
  clust_net <- reactive({
    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees
    medK <- as.numeric(sn$sub_nets$median)

    clust_net <- list() 
    clust_net[["genes"]] <- cluster_coexp( sub_net$genes, medK = medK, flag_plot = FALSE )
    return(clust_net)

  })
  
  
  # Output of subnetowrk table
  observeEvent(
    input$generate_subnet, 
    {output$subnetwork <- renderTable(sn$sub_nets)}
  )
  
  







  # CLUSTER GENES

  observeEvent(
    {input$run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees
      medK <- as.numeric(sn$sub_nets$median)
     

      # network output
      output$CNtext = renderText(
        # input$run
        # req(input$run) # to prevent print at first lauch
        # isolate(h6("Network of Clustered Genes"))
        "Network of Clustered Genes"
      )
      output$network <- renderPlot(
        {plot_network( sub_net$gene, clust_net()$gene, medK )},
        width = 500,
        height = 500
      )


      # heatmap output
      output$CHtext = renderText({
        input$run
        req(input$run) # to prevent print at first lauch
        isolate(h6("Heatmap of Clustered Genes"))
      })
      output$heatmap <- renderPlot(
        {plot_coexpression_heatmap( sub_net$gene, clust_net()$gene, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )


      # binarized heatmap output
      output$CHBtext = renderText({
        input$run
        req(input$run) # to prevent print at first lauch
        isolate(print("Binarized Heatmap of Clustered Genes"))
      })
      output$Bheatmap <- renderPlot(
        {plot_coexpression_heatmap( sub_net$gene, clust_net()$gene )},
        width = 500,
        height = 500
      )
      

      # clustering genes table output
      output$CG_table_text = renderText({
        # input$runGC
        # req(input$runGC) # to prevent print at first lauch
        # isolate(print("Filtered Genes"))
        "Clustering Genes"
      })
      output$CG_table <- renderDT(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


    }
  )



  # GENE CONNECTIVITY

  observeEvent(
    {input$runGC},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)
      m <- match(clust_net()$genes$clusters$genes, rownames(sub_net$genes))


      # density output
      output$GCdensityGtext = renderText({
        input$runGC
        req(input$runGC) # to prevent print at first lauch
        isolate(print("Density plot of Gene Connectivity"))
      })
      output$GCdensityG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )


      # histogram output
      output$GChistogramGtext = renderText({
        input$runGC
        req(input$runGC) # to prevent print at first lauch
        isolate(print("Histogram of Gene Connectivity"))
      }) 
      output$GChistogramG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xybreaks = input$xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )


      # density output - subset by clusters
      output$GCdensitySubsetGtext = renderText({
        input$runGC
        req(input$runGC) #to prevent print at first lauch
        isolate(print("Density plot of Gene Connectivity subset by their clusters"))
      })
      output$GCdensitySubsetG <- renderPlot(
        { plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "density")},
        width = 500,
        height = 500
      )


      # histogram output - subset by clusters
      output$GChistogramSubsetGtext = renderText({
        input$runGC
        req(input$runGC) #to prevent print at first lauch
        isolate(print("Histogram of Gene Connectivity subset by their clusters"))
      }) 
      output$GChistogramSubsetG <- renderPlot(
        { plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xybreaks = input$xybreaks,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "hist")},
        width = 500,
        height = 500
      )


    }
  )


  # FUNCTIONAL OUTLIERS

  observeEvent(
    {input$runFO},
    {
      sub_net <- sn$sub_nets$sub_net
      filt_min <- 6 
      clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))


      # unselected genes table output
      output$FO_unselected_text = renderText({
        # input$runGC
        # req(input$runGC) # to prevent print at first lauch
        # isolate(print("Filtered Genes"))
        "Unselected Genes"
      })
      output$genes_unselected_table <- renderDT(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


      # selected genes table output
      output$FO_selected_text = renderText({
        # input$runGC
        # req(input$runGC) # to prevent print at first lauch
        # isolate(print("Unfiltered Genes"))
        "Selected Genes"
      })
      output$genes_selected_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
      
    }
  )

  # ERROR MESSAGES

  output$CG_error <- renderText({
    print("Please upload/generate a gene list in OPTIONS")
  })

  output$GC_error <- renderText({
    print("Please upload/generate a gene list in OPTIONS")
  }) 

  output$FO_error <- renderText({
    print("Please upload/generate a gene list in OPTIONS")
  }) 

}


