ui <- fluidPage(
  titlePanel(title=div(img(src="ODClogo.png", height = 50), "OutDeCo")),
  navbarPage("",
             tabPanel(
               title="Home",
               
               navlistPanel(
                 tabPanel(
                   title="About",
                   "What is OutDeCo?"       
                 ),
                 
                 tabPanel(
                   title="Differential Analysis",
                   "information about DE",
                 ),
                 
                 
                 tabPanel(
                   title="Cluster Genes",
                   "information about Cluster genes",
                 ),
                 
                 tabPanel(
                   title="Gene Connectivity",
                   "information about Gene Connectivity",
                 ),
                 
                 
                 tabPanel(
                   title="Functional Outliers",
                   "information about Functional Outliers",
                 ),
                 
                 tabPanel(
                   title="Gene Set Enrichment Analysis",
                   "information about Enrichment Analysis",
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
               #Sidebar
               navlistPanel(
                 tabPanel(
                   title="Cluster Genes",
                   "Cluster genes Page",
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
