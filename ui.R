library(shinyWidgets)
library(shinythemes)
library(bslib)
library(DT)


ui <- fluidPage(
 
  titlePanel(title=div(img(src="ODClogo.png", height = 80), "OutDeCo")),
  theme = bs_theme(version = 5, bootswatch = "sandstone"),
  
  navbarPage("", collapsible = TRUE,
          
             tabPanel(
               title="Home",
               icon = icon("home"),

               navlistPanel(
                 id = "Header", selected = NULL, well = FALSE, fluid = FALSE, widths = c(3, 8) ,
             
                 tabPanel(title="About",
                   h3("What is OutDeCo?"),
                   h6("Outlier detection through co-expression. Using the tools on this website you can:"),
                   p("- Run a differential expression (DE) analysis "),
                   p("- Assess your DE results using gene co-expression properties"),
                   p("- Report a functional outlier assessment"),
                   p("- Run a network connectivity analysis of DE results within a gene co-expression network"),
                 ),
                 
                 tabPanel(title="Differential Analysis",
                   h3("Differential Analysis Methods"),
                   h5("Methods:"),
                   h6("wilcox"),
                   p("Compares means of two groups to analyse count data and test for differential expression."),
                   h6("DESeq"),
                   p("Uses geometric normalisation to analyse count data and test for differential expression."),
                   h6("edgeR"),
                   p("Uses weighted mean of log ratio to analyse count data and test for differential expression."),
                   br(),
                   em("Note: This tool is under construction")
                   
                  ),
                 
                 
                 tabPanel(title="Cluster Genes",
                   h3("Cluster Genes"),
                   p('Creates modules which are clusters of genes that are hightly co-expressed'),
                   h5("Plot Types"),
                   h6("Heatmap"),  
                   p("up-regulated or down-regulated network plot"),
                   h6("Network"),  
                   p("up regulated or down-regulated network plot where nodes are genes and the weight of the edges corresponds to the co-expression levels between its endpoints"),
                   h6("Binarized heatmap"),  
                   p("up-regulated or down-regulated binary co-expression sub-network"),
                   br(),
                   em("Note: This tool is under construction")
                 ),
                 
                 tabPanel(title="Gene Connectivity",
                   h3("Cluster Genes"),
                   p('Creates modules which are clusters of genes that are hightly co-expressed'),
                   h5("Plot Types"),
                   h6("Heatmap"),  
                   p("up-regulated or down-regulated network plot"),
                   h6("Network"),  
                   p("up regulated or down-regulated network plot where nodes are genes and the weight of the edges corresponds to the co-expression levels between its endpoints"),
                   h6("Binarized heatmap"),  
                   p("up-regulated or down-regulated binary co-expression sub-network"),
                   
                   br(),
                   em("Note: This tool is under construction"),
                   
                 ),
                 
                 
                 tabPanel(
                   title="Functional Outliers",
                   "Functional outliers are genes that have been identified to be potentially dysregulated. 
                   They are the genes that are Differentially Expressed but do not show local co-expression",
                 
                   h5("Analysis Types"),
                   h6("Coexpression heatmap"),
                   "- Up and Down regulated",
                   h6("Network"),
                   "- Up and Down regulated",
                   
                   br(),
                   em("Note: This tool is under construction"),
                   
                 ),
                 
                 tabPanel(
                   title="Gene Set Enrichment Analysis",
                   "information about Enrichment Analysis",
      
                   br(),
                   em("Note: This tool is under construction")
                 ),
                 
               )
             ),
             
             
             tabPanel(
               title="Run DE",
               navlistPanel(
                 tabPanel(
                   title="wilcox",
                   "wilcox placeholder",
                 ),
                 
                 tabPanel(
                   title="DESeq",
                   "DESeq placeholder",
                 ),
                 
                 
                 tabPanel(
                   title="edgeR",
                   "edgeR Placeholder",
                 ),
               ),
               
               
             ),
             tabPanel(
               title="Assess DE",
               dropdown(
                 
                 tags$h3("List of Input"),
                 
                 pickerInput(inputId = 'xcol2',
                             label = 'X Variable',
                             choices = names(iris),
                             options = list(`style` = "btn-info")),
                 style = "gradient", icon = icon("cog"),
                 status = "primary", width = "300px",
                 animate = animateOptions(
                   enter = animations$fading_entrances$fadeInLeftBig,
                   exit = animations$fading_exits$fadeOutLeftBig
                 )
               ),
               
              #Sidebar
              #DROPDOWN
              
              navlistPanel(
                tabPanel(
                  title="Cluster Genes",
                  "Cluster genes Page",
                  #Navigation Bar for types of plots
                  tabsetPanel(
                    tabPanel(
                      title="Plot 1", 
                    ),
                    tabPanel(
                      title ="Plot 2"
                    ),
                    tabPanel(
                      title="Plot 3"
                    )
                  ),
                ),
                
                tabPanel(
                  title="Gene Connectivity",
                  "Gene Connectivity Page",
                ),
                 
                 
                 tabPanel(
                   title="Functional Outliers",
                   "Functional Outliers Page",
                 ),
                 
                 tabPanel(
                   title="Gene Set Enrichment Analysis",
                   "GSE Page",
                 ),
                 
               ),
               
             )
             
             
  ),
)
