source("./src/plot_functions.R", local = TRUE)
source("./src/cluster_coexp.R", local = TRUE)
source("./src/subset_network_hdf5.R", local = TRUE)
source("./src/calc_DE.R", local = TRUE)

# Warnings silenced for wilcox
options(warn=-1)
defaultW <- getOption("warn")


server <- function(input, output, session) {

  # Removing elements that are not functional without subnetwork
  hide(id = "GC_dropdown")
  hide(id = "GC_dropdown_DE")
  hide(id = "CG_dropdown")
  hide(id = "CG_dropdown_DE")
  hide(id = "FO_dropdown")
  hide(id = "FO_dropdown_DE")
  hide(id = "assess_run_de")


  hide(id = "sepLabelsButton")
  hide(id = "sepCountsButton")

  #Run Differential Expression
  hide(id="DE_options")
  hide(id="run_DE")



  # PLOTS/TABLES HEADERS
  # Run DE
  hide(id = "vol")
  hide(id = "MA")
  # Clustering 
  hide(id="CG_network_text")
  hide(id="CG_heatmap_text")
  hide(id="CG_bheatmap_text")
  hide(id="CG_table_text")
  hide(id="CGupreg_network_text")
  hide(id="CGupreg_heatmap_text")
  hide(id="CGupreg_bheatmap_text")
  hide(id="CGdownreg_network_text")
  hide(id="CGdownreg_heatmap_text")
  hide(id="CGdownreg_bheatmap_text")


  # Gene Connectivity
  hide(id="GCdensityG_text")
  hide(id="GChistogramG_text")
  hide(id="GCdensitySubsetG_text")
  hide(id="GChistogramSubsetG_text")
  hide(id="GCdensityG_upreg_text")
  hide(id="GChistogramG_upreg_text")
  hide(id="GCdensitySubsetG_upreg_text")
  hide(id="GChistogramSubsetG_upreg_text")
  hide(id="GCdensityG_downreg_text")
  hide(id="GChistogramG_downreg_text")
  hide(id="GCdensitySubsetG_downreg_text")
  hide(id="GChistogramSubsetG_downreg_text")

  #Functional Outliers
  hide(id="FO_network_text")
  hide(id="FO_heatmap_text")
  hide(id="FOnetwork_upreg_text")
  hide(id="FOheatmap_upreg_text")
  hide(id="FOnetwork_downreg_text")
  hide(id="FOheatmap_downreg_text")
  hide(id="genes_not_keep_table_text")
  hide(id="genes_keep_table_text")





  ##########################################################################################
  #                                                                                        #
  #                                    RUN DE                                              #
  #                                                                                        #
  ##########################################################################################
  
  
  #####################  UPLOAD COUNTS DATA ###########################

  # Make countsData
  countsData <- reactive({
    ServerCountsFile <- input$counts_file
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    req(ServerCountsFile)
    validate(need(extCountsFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extCountsFile)) {
      return ()
    }
    if (extCountsFile == "csv") {
      read.table(file=ServerCountsFile$datapath, sep=input$sepCountsButton, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerCountsFile$datapath, sep=input$sepCountsButton, header=TRUE) 
    }
    
    
  })

  observe({
     # counts_file from fileInput() function
    ServerCountsFile <- req(input$counts_file)
    
  #   # extensions tool for format validation
    extCountsFile <- tools::file_ext(ServerCountsFile$datapath)
    if (is.null(input$counts_file)) {
      return ()
    } else{
      if (extCountsFile == "txt") {
        label = paste("Delimiters for", extCountsFile, "file")
        choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
      } else if (extCountsFile == "tsv") {
        label = paste("Delimiter: Tab")
        choice <- (Tab="\t")
      } else {
        label = paste("Delimiter: Comma")
        choice <- (Comma=",")
      }
      updateRadioButtons(session, "sepCountsButton", label = label, choices = choice)
      }

      #print(counts_data[:])
    })

    output$UICountsContent <- renderDataTable(
        countsData(), options = list(
          pageLength = 25
        )
      )
      
    observeEvent(input$counts_file, {
      show(id = "sepCountsButton")
    })

  ########################### UPLOAD LABELS DATA ###########################
  # Make labelsData
  labelsData <- reactive({
    ServerLabelsFile <- input$labels_file
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    req(ServerLabelsFile)
    validate(need(extLabelsFile == c("csv", "tsv", "txt"), "Please upload a csv, tsv or txt file."))
    if (is.null(extLabelsFile)) {
      return ()
    }
    if (extLabelsFile == "csv") {
      read.table(file=ServerLabelsFile$datapath, sep=input$sepLabelsButton, header=TRUE, row.names = 1)
    } else {
      read.table(file=ServerLabelsFile$datapath, sep=input$sepLabelsButton, header=TRUE)
    }
    
    
  })

  observeEvent(input$labels_file, {
    show(id = "sepLabelsButton")
  })



  observe({
     # labels_file from fileInput() function
    ServerLabelsFile <- req(input$labels_file)
    
  #   # extensions tool for format validation
    extLabelsFile <- tools::file_ext(ServerLabelsFile$datapath)
    if (is.null(input$labels_file)) {
      return ()
    } else {
      if (extLabelsFile == "txt") {
        label = paste("Delimiters for", extLabelsFile, "file")
        choice <-c(Comma=",", Semicolon=";", Tab="\t", Space=" ")
      } else if (extLabelsFile == "tsv") {
        label = paste("Delimiter: Tab")
        choice <- (Tab="\t")
      } else {
        label = paste("Delimiter: Comma")
        choice <- (Comma=",")
      }
      updateRadioButtons(session, "sepLabelsButton", label = label, choices = choice)
      }
    })


  # handles rendering DT table of labels file

  output$UILabelContent <- renderDataTable(
    labelsData(), options = list(
      pageLength = 100
    )
  )

  # rendering DT table for RUN DE options (to select cases)
  output$UILabelContentSelection <- renderDataTable(
    labelsData(), options = list(
      pageLength = 100
    )
  )

    # rendering DT table for RUN DE options (to remove cases)
  output$UILabelContentRemoveSelection <- renderDataTable(
    labelsData(), options = list(
      pageLength = 100
    ) 
  )

  # Plot the data
  observeEvent(input$case_control_method, {
      options <- names(labelsData())
      updateSelectInput(session, 
        inputId="select_column",
        "Select column to group", 
        choices = options[2:length(options)], 
        selected = NULL
      )
      show(id="run_DE")
  })

  # Group by label option 
  observeEvent(input$select_column, {
      # Update the case selection with levels of selected column 
      var <- labelsData()[[input$select_column]]
      lvl <- levels(as.factor(var))
      updateSelectInput(session, 
        inputId="select_case", 
        "Select case to analyse", 
        choices = lvl, 
        selected = NULL
      )
  })

  # countsData <- reactive({
  #   if ( is.null(input$counts_file)) return(NULL)
  #   inFile <- input$counts_file
  #   file <- inFile$datapath
  #   # load the file into new environment and get it from there
  #   e = new.env()
  #   name <- load(file, envir = e)
  #   data <- e[[name]]
  # })

  case_selected <- reactive({
    input$UILabelContentSelection_rows_selected
  })

  conditions_selected <- reactive({
    input$UILabelContentRemoveSelection_rows_selected
  })

  
    
  # Switch to labels tab if labels file is uploaded
  observeEvent(input$labels_file, {
    updateTabsetPanel(session, "counts_labels_tabset", selected = "Labels File")
  })
  # Switch to counts tab if counts file is uploaded
  observeEvent(input$counts_file, {
    updateTabsetPanel(session, "counts_labels_tabset", selected = "Counts File")
  })




  ########################### RUN DE ###########################
  observe({
    if (!is.null(input$labels_file) && !is.null(input$labels_file)) {
      show(id="DE_options")
      hide(id="runDE_error")
    }
  })
  
  de <- reactiveValues(
    deg_output = NULL, 
  )


  observeEvent(input$run_DE, {
    labels <- labelsData()
    counts_data <- countsData()
    deg <- NULL

    # var <- labelsData()[[input$select_column]]
    if (input$case_control_method == "Choose Case by Label") {
      var <- input$select_column
      case <- input$select_case

      # Format labels$var
      labels_var <- labels[[paste0(var)]]

      #Initialise the variables of the chosen column to all be 1
      groups <- rep(1, length(labels_var))
      
      # Pick the case, relabel as 2
      groups[labels_var == case] = 2   


      filt = groups != 0 
      deg <- calc_DE(counts_data[,filt], groups[filt], input$DE_method) 
      de$deg_output <- deg

    } else {
      cases <- case_selected()
      conditions <- conditions_selected()
      
      #Initalise all values to 1
      groups <- rep(0, nrow(labels))

      for (c in cases) {
        groups[c] = 2
      }
      
      for (d in conditions) {
        groups[d] = 1
      }
      # No cases have been selected
      filt = groups != 0 
      if (is.null(cases)) {
        shinyalert(title = "Invalid Input", text = "Please select cases to assess", type = "error")
      } else if (is.null(conditions)) {
        shinyalert(title = "Invalid Input", text = "Please select conditions to assess", type = "error")
      # Valid input - cases and control selected
      } else {
        deg <- calc_DE(counts_data[,filt], groups[filt], input$DE_method) 
        de$deg_output <- deg
      }
      

    }
    

    # Volcano Plot
    if (!is.null(deg)) {
      show(id="vol")
      output$DEplot <- renderPlot(
              {plot( deg$degs$log2_fc, -log10(deg$degs$pvals),  
              pch=19, bty="n", 
              xlab="log2 FC", ylab="-log10 p-vals" )},
              width = 450,
              height = 450
      )

      #MA Plot
      show(id="MA")
      #output$DE_MA_text = renderText("MA Plot")
      output$DEplot_average <- renderPlot(
              {plot( log2(deg$degs$mean_cpm),  deg$degs$log2_fc,  
              pch=19, bty="n", 
              ylab="log2 FC", xlab="Average expression (log2 CPM + 1)")},
              width = 450,
              height = 450
      )
    show(id = "assess_run_de")
    }
    
    }
  )

  observeEvent(input$assess_run_de, { 
    updateTabsetPanel(session, inputId="navpage", selected="Assess DE")
    updateTabsetPanel(session, "subnetwork_file_tabset_DE", selected = "Subnetwork")
    sn$sub_nets <- NULL
    output$DE_table <- renderDataTable(
        {de$deg_output$degs},
    )
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
  
  ##########################################################################################
  #                                                                                        #
  #                                    ASSESS DE DATA                                      #
  #                                                                                        #
  ##########################################################################################

  # sub_nets
  sn <- reactiveValues(
    sub_nets_DE = NULL,
    sub_nets = NULL, 
  )

  observeEvent(input$generate_subnet_DE, {
    if (is.null(de$deg_output)) {
      shinyalert(title = "Invalid Input", text = "Please first Run DE", type = "error")
      updateTabsetPanel(session, inputId="navpage", selected="Run DE")
    } else {
      sn$sub_nets_DE <- subset_network_hdf5(de$deg_output$degs, tolower(input$network_type), dir="../networks/")
      show(id = "CG_dropdown_DE")
      hide(id = "CG_error_DE")
      show(id = "GC_dropdown_DE")
      hide(id = "GC_error_DE")
      show(id = "FO_dropdown_DE")
      hide(id = "FO_error_DE")
    }
  })

  observeEvent(
    input$generate_subnet_DE, 
    {output$subnetwork_DE <- renderTable(sn$sub_nets_DE)}
  )

  clust_net_DE <- reactive({

    sub_net <- sn$sub_nets_DE$sub_net
    node_degrees <- sn$sub_nets_DE$node_degrees
    medK <- as.numeric(sn$sub_nets_DE$median)

    clust_net_DE <- list() 
    
    # For DE data 
    deg_sig <- sn$sub_nets_DE$deg_sig
    fc_sig  <- sn$sub_nets_DE$fc_sig
    clust_net_DE[["down"]]  <- cluster_coexp(sub_net$down, medK = medK, flag_plot = FALSE)
    clust_net_DE[["up"]]  <- cluster_coexp( sub_net$up, medK = medK, flag_plot = FALSE)

    return(clust_net_DE)
  })


  ################################ CLUSTER GENES ########################################

  observeEvent(
    {input$runCGDE},
    {
      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees
      medK <- as.numeric(sn$sub_nets_DE$median)

      
      # upregulated network 
      show(id="CGupreg_network_text")
      output$upregNetwork <- renderPlot(
        {plot_network(sub_net$up, clust_net_DE()$up, medK)}, 
        width = 500, 
        height = 500 
      )

      # upregulated heatmap 
      show(id="CGupreg_heatmap_text")
      output$upregHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up, flag_plot_bin = FALSE)}, 
        width = 500,
        height = 500
      )

      # upregulated binarized heatmap 
      show(id="CGupreg_bheatmap_text")
      output$upregbinHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up)}, 
        width = 500, 
        height = 500
      )
      
      # downregulated network 
      show(id="CGdownreg_network_text")
      output$downregNetwork <- renderPlot(
        {plot_network(sub_net$down, clust_net_DE()$down, medK)},
        width = 500, 
        height = 500
      )

      # downregulated heatmap
      show(id="CGdownreg_heatmap_text")
      output$downregHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down, flag_plot_bin = FALSE)}, 
        width = 500, 
        height = 500 
      )

      # downregulated binarized heatmap
      show(id="CGdownreg_bheatmap_text")
      output$downregbinHeatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down)}, 
        width = 500, 
        height = 500
      )

      # clustering genes table output
      # show(id="CG_table_text")
      # output$CG_table <- renderDataTable(
      #   {EGAD::attr.human[match(clust_net()$genes$clusters$genes,EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome],input$chooseGeneNo),]},
      #   # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      # )


    }
  )

  ################################ GENE CONNECTIVITY ######################################

  observeEvent(
    {input$runGCDE},
    {
      
      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees  
      medK <- as.numeric(sn$sub_nets_DE$median)
      

      # density - upreg
      show(id="GCdensityG_upreg_text")
      output$GCdensityGupreg <- renderPlot(
        {plot_scatter(node_degrees$up[,1]/node_degrees$n_genes_total, 
                    node_degrees$up[,2]/node_degrees$n_genes_up, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )

      # histogram - upreg 
      show(id="GChistogramG_upreg_text")
      output$GChistogramGupreg <- renderPlot(
        {plot_scatter(node_degrees$up[,1]/node_degrees$n_genes_total, 
                    node_degrees$up[,2]/node_degrees$n_genes_up, 
                    xybreaks = input$xybreaks_DE,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      # density by subsets - upreg 
      show(id="GCdensitySubsetG_upreg_text")
      output$GCdensitySubsetGupreg <- renderPlot(
        {m <- match(clust_net_DE()$up$clusters$genes, rownames(sub_net$up))
         plot_scatter(node_degrees$up[m,1]/node_degrees$n_genes_total, 
                      node_degrees$up[m,2]/node_degrees$n_genes_up, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$up$clusters, flag = "density")},
        width = 500,
        height = 500
      )

      # histogram by subsets - upreg 
      show(id="GChistogramSubsetG_upreg_text")
      output$GChistogramSubsetGupreg <- renderPlot(
        {m <- match(clust_net_DE()$up$clusters$genes, rownames(sub_net$up))
         plot_scatter(node_degrees$up[m,1]/node_degrees$n_genes_total, 
                      node_degrees$up[m,2]/node_degrees$n_genes_up, 
                      xybreaks = input$xybreaks_DE,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$up$clusters, flag = "hist")},
        width = 500,
        height = 500
      )

      # density - downreg
      show(id="GCdensityG_downreg_text")
      output$GCdensityGdownreg <- renderPlot(
        {plot_scatter(node_degrees$down[,1]/node_degrees$n_genes_total, 
                    node_degrees$down[,2]/node_degrees$n_genes_down, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )

      # histogram - downreg 
      show(id="GChistogramG_downreg_text")
      output$GChistogramGdownreg <- renderPlot(
        {plot_scatter(node_degrees$down[,1]/node_degrees$n_genes_total, 
                    node_degrees$down[,2]/node_degrees$n_genes_down, 
                    xybreaks = input$xybreaks_DE,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      # density by subsets - downreg 
      show(id="GCdensitySubsetG_downreg_text")
      output$GCdensitySubsetGdownreg <- renderPlot(
        {m <- match(clust_net_DE()$down$clusters$genes, rownames(sub_net$down))
         plot_scatter(node_degrees$down[m,1]/node_degrees$n_genes_total, 
                      node_degrees$down[m,2]/node_degrees$n_genes_down, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$down$clusters, flag = "density")},
        width = 500,
        height = 500
      )

      # histogram by subsets - downreg 
      show(id="GChistogramSubsetG_downreg_text")
      output$GChistogramSubsetGdownreg <- renderPlot(
        {m <- match(clust_net_DE()$down$clusters$genes, rownames(sub_net$down))
         plot_scatter(node_degrees$down[m,1]/node_degrees$n_genes_total, 
                      node_degrees$down[m,2]/node_degrees$n_genes_down, 
                      xybreaks = input$xybreaks_DE,
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net_DE()$down$clusters, flag = "hist")},
        width = 500,
        height = 500
      )

    }
  )

  ################################ FUNCTIONAL OUTLIERS ######################################
  observeEvent(
    {input$runFODE},
    {

      sub_net <- sn$sub_nets_DE$sub_net
      node_degrees <- sn$sub_nets_DE$node_degrees  
      medK <- as.numeric(sn$sub_nets_DE$median)

      filt_min <- input$filtmin

      show(id="FOheatmap_upreg_text")
      output$FOheatmap_upreg <- renderPlot(
        {plot_coexpression_heatmap(sub_net$up, clust_net_DE()$up, filt = TRUE, flag_plot_bin = FALSE)}, 
        width = 500,
        height = 500 
      )

      show(id="FOnetwork_upreg_text")
      output$FOnetwork_upreg <- renderPlot(
        {plot_network(1-sub_net$up, clust_net_DE()$up, 1 - medK)}, 
        width = 500, 
        height = 500
      )

      show(id="FOheatmap_downreg_text")
      output$FOheatmap_downreg <- renderPlot(
        {plot_coexpression_heatmap(sub_net$down, clust_net_DE()$down, filt = TRUE, flag_plot_bin = FALSE)}, 
        width = 500, 
        height = 500 
      )

      show(id="FOnetwork_downreg_text")
      output$FOnetwork_downreg <- renderPlot(
        {plot_network(1 - sub_net$down, clust_net_DE()$down, 1 - medK)}, 
        width = 500, 
        height = 500
      )

      # # genes in module table output
      # show(id="genes_not_keep_table_text")
      # output$genes_not_keep_table <- renderDataTable(
      #   { clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      #     clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      #     genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
      #     EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
      #   # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      # )


      # # functional outliers table output
      # show(id="genes_keep_table_text")
      # output$genes_keep_table <- renderDataTable(
      #   { clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      #     clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      #     genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
      #     EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
      #   # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      # )
      
    }
  )
  

  ##########################################################################################
  #                                                                                        #
  #                                  ASSESS GENE LIST                                      #
  #                                                                                        #
  ##########################################################################################

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

    read.table(file=ServerDEFile$datapath, sep=input$sepButton, header=TRUE)
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

  observeEvent(input$gene_list_selection, {
    if (input$gene_list_selection == "Upload Gene List") {
      updateTabsetPanel(session, "subnetwork_file_tabset", selected = "File")
    } else {
      updateTabsetPanel(session, "subnetwork_file_tabset", selected = "Subnetwork")
    }
  })

  observeEvent(input$generate_subnet, {
    gene_list <- NULL
    if (is.null(input$gene_list_selection)) {
      shinyalert(title = "Invalid Input", text = "Please choose a gene list method", type = "error")
    } else {
      # GENERATE GENE LIST
      if (input$gene_list_selection == "Generate Gene List") {

        if (str_detect(input$chooseChrome, "chr[XY]") == FALSE && str_detect(input$chooseChrome, "chr[0-9]") == FALSE) {
          shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
          gene_list <- NULL

        } else if (str_detect(substring(input$chooseChrome,4), "[0-9]")) {
          if (strtoi(substring(input$chooseChrome,4)) < 1 || strtoi(substring(input$chooseChrome,4)) > 22) {
            shinyalert(title = "Invalid Input", text = "Please enter a Chromosome between 1 - 22, X, Y", type = "error")
            gene_list <- NULL

          } else {
            gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo,)
            print(gene_list)
          }

        } else if (input$chooseGeneNo == "" || input$chooseGeneNo < 0) { 
          shinyalert(title = "Invalid Input", text = "Please enter a valid number of Genes", type = "error")
          gene_list <- NULL

        } else { 
          gene_list <- sample( EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo,)

        }
          
      } else {
        # Invalid Input - user hasn't uploaded file
        if (is.null(input$DEFile)) {
          shinyalert(title = "Invalid Input", text = "Please upload a gene list file", type = "error")
        } else {
          gene_list <- read.delim(file = input$DEFile$datapath, header = FALSE, sep = "\n", dec = ".")[,1]
        }

      }
      #Valid Input
      if (!is.null(gene_list)) { 
        sn$sub_nets <- subset_network_hdf5_gene_list(gene_list, tolower(input$network_type), dir="../networks/")
        show(id = "CG_dropdown")
        hide(id = "CG_error")
        show(id = "GC_dropdown")
        hide(id = "GC_error")
        show(id = "FO_dropdown")
        hide(id = "FO_error")
        # Clear data
        output$network <- NULL
        output$heatmap <- NULL
        output$Bheatmap <- NULL
        output$CG_table <- NULL
        output$GCdensityG <- NULL
        output$GChistogramG <- NULL
        output$GCdensitySubsetG<- NULL
        output$GChistogramSubsetG <- NULL
        output$FO_heatmap <- NULL
        output$FO_network <- NULL
        output$genes_not_keep_table <- NULL
        output$genes_keep_table <- NULL
        # Reset Checkboxes
        updateAwesomeCheckboxGroup(
          inputId = "clusterPlotOptions_genelist",
          choices = c("Network", "Heatmap", "Binarized Heatmap"),
          status = ""
        )
        updateAwesomeCheckboxGroup(
          inputId = "GCPlotOptions_genelist",
          choices = c("Density", "Histogram", "Clustered Density", "Clustered Histogram"),
          status = ""
        )
        updateAwesomeCheckboxGroup(
          inputId = "FOPlotOptions_genelist",
          choices = c("Network", "Heatmap"),
          status = ""
        )
        updateAwesomeCheckboxGroup(
          inputId = "FO_table_options",
          choices = c("Functional Outliers", "Genes in Module"),
          status = ""
        )
      } 
    }
  })

  observeEvent(
    input$generate_subnet, 
    {output$subnetwork <- renderTable(sn$sub_nets)}
  )

  clust_net <- reactive({

    sub_net <- sn$sub_nets$sub_net
    node_degrees <- sn$sub_nets$node_degrees
    medK <- as.numeric(sn$sub_nets$median)

    clust_net <- list() 
    clust_net[["genes"]] <- cluster_coexp(sub_net$genes, medK = medK, flag_plot = FALSE)
    
    return(clust_net)

  })

  ##################### CLUSTER GENES #####################

  observeEvent(
    {input$run},
    {
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees
      medK <- as.numeric(sn$sub_nets$median)

      # network output
      show(id="CG_network_text")
      CG_network <- function(){plot_network(sub_net$genes, clust_net()$genes, medK)}
      output$network <- renderPlot(
        {CG_network()},
        width = 500,
        height = 500
      )


      # heatmap output
      show(id="CG_heatmap_text")
      CG_heatmap <- function(){plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, flag_plot_bin = FALSE)}
      output$heatmap <- renderPlot(
        {CG_heatmap()},
        width = 500,
        height = 500
      )


      # binarized heatmap output
      show(id="CG_bheatmap_text")
      output$Bheatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes)},
        width = 500,
        height = 500
      )

      # clustering genes table output
      show(id="CG_table_text")
      output$CG_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes,EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome],input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )

      output$downloadCG_genelist <- downloadHandler(
        contentType = "image/png",
        filename = function() {
          paste("ClusteredNetwork", ".png")
        },
        content = function(file) {
          png(file, width=500, height=500)
          CG_network()
          dev.off()
        }
      )

      output$downloadCG_genelist <- downloadHandler(
        contentType = "image/png",
        filename = function() {
          paste("ClusteredHeatmap", ".png")
        },
        content = function(file) {
          png(file, width=500, height=500)
          CG_heatmap()
          dev.off()
        }
      )
    }
  )

  #Download Plots
  


 
   

    
  



  
  #Download Tables
    # filename = function() {
    #   paste("data-", Sys.Date(), ".csv", sep="")
    # },
    # content = function(file) {
    #   write.csv(data, file)
    # }


  ##################### GENE CONNECTIVITY #####################

  observeEvent(
    {input$runGC},
    {
      
      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)
      m <- match(clust_net()$genes$clusters$genes, rownames(sub_net$genes))

      # density output
      show(id="GCdensityG_text")
      output$GCdensityG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "density")},
        width = 500,
        height = 500
      )


      # histogram output
      show(id="GChistogramG_text")
      output$GChistogramG <- renderPlot(
        {plot_scatter(node_degrees$genes[,1]/node_degrees$n_genes_total, 
                    node_degrees$genes[,2]/node_degrees$n_genes, 
                    xybreaks = input$xybreaks,
                    xlab="Global node degree", 
                    ylab="Local node degree", flag= "hist")},
        width = 500,
        height = 500
      )

      show(id="GCdensitySubsetG_text")
      # density output - subset by clusters
      output$GCdensitySubsetG <- renderPlot(
        {plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
                      node_degrees$genes[m,2]/node_degrees$n_genes, 
                      xlab="Global node degree", 
                      ylab="Local node degree", 
                      clusters = clust_net()$genes$clusters, flag = "density")},
        width = 500,
        height = 500
      )


      # histogram output - subset by clusters
      show(id="GChistogramSubsetG_text")
      output$GChistogramSubsetG <- renderPlot(
        {plot_scatter(node_degrees$genes[m,1]/node_degrees$n_genes_total, 
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


  ##################### FUNCTIONAL OUTLIERS #####################

  observeEvent(
    {input$runFO},
    {

      sub_net <- sn$sub_nets$sub_net
      node_degrees <- sn$sub_nets$node_degrees  
      medK <- as.numeric(sn$sub_nets$median)

      filt_min <- input$filtmin

      clust_size <- plyr::count(clust_net()$genes$clusters$labels)
      clust_keep <- clust_size[clust_size[,2] < filt_min ,1]
      genes_keep <- !is.na(match(clust_net()$genes$clusters$labels, clust_keep))
     

      # heatmap output
      show(id="FO_heatmap_text")
      output$FO_heatmap <- renderPlot(
        {plot_coexpression_heatmap(sub_net$genes, clust_net()$genes, filt = TRUE, flag_plot_bin = FALSE)},
        width = 500,
        height = 500
      )

      # network output
      show(id="FO_network_text")
      output$FO_network <- renderPlot(
        {plot_network(1-sub_net$genes, clust_net()$genes, 1 - medK)},
        width = 500,
        height = 500
      )

      # genes in module table output
      show(id="genes_not_keep_table_text")
      output$genes_not_keep_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[!genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )


      # functional outliers table output
      show(id="genes_keep_table_text")
      output$genes_keep_table <- renderDataTable(
        {EGAD::attr.human[match(clust_net()$genes$clusters$genes[genes_keep],EGAD::attr.human$name[EGAD::attr.human$chr==input$chooseChrome], input$chooseGeneNo),]},
        # options=list(columnDefs = list(list(visible=FALSE, targets=c(0,1,2,3))))
      )
      
    }
  )




  ##################### ERROR MESSAGES #####################

  output$CG_error <- renderText({
    print("Please upload/generate a gene list in NETWORK OPTIONS")
  })

  output$CG_error_DE <- renderText({
    print("Please upload/use DE Data in NETWORK OPTIONS")
  })

  output$GC_error <- renderText({
    print("Please upload/generate a gene list in NETWORK OPTIONS")
  }) 

  output$GC_error_DE <- renderText({
    print("Please upload/use DE Data in NETWORK OPTIONS")
  })

  output$FO_error <- renderText({
    print("Please upload/generate a gene list in NETWORK OPTIONS")
  }) 

  output$FO_error_DE <- renderText({
    print("Please upload/use DE Data in NETWORK OPTIONS")
  })

}


