# Load packages
library(shiny)
library(shinythemes)
library(mgcv)
library(shinycssloaders)
library(dplyr)

##################################################
## Define UI for app
ui <- fluidPage(
  theme = shinytheme("flatly"),
  title = "PANORAMA",
  
  br(),
  
  ##################################################
  ## Sidebar layout
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Title
      div(p("PANORAMA"),
          align = "center",
          style = "font-size:32px;font-weight:bold;"),
      
      br(),
      
      # View type
      selectInput("plotType", "View type",
                c("Pan-cancer" = "pancancer",
                "Single model" = "singleTT")),
      
      # Gene name
      textInput("gene", "Gene name", 
                value = "ERBB2"),
      
      # Sort order
      conditionalPanel(condition = "input.plotType == 'pancancer'",
        selectInput("sortOrder", "Sort",
                  c("Individually" = "individually",
                    "CNA" = "bycna",
                    "Methylation" = "bymeth",
                    "Combined" = "bycombined",
                    "Alphabetical" = "alphabetical"))),
              
      # Show combined?
      conditionalPanel(condition = "input.plotType == 'pancancer'",
        div(checkboxInput(inputId = "showCombined",
                          label = "Show combined model️︎"), 
            align = "right")),
      
      # Tumor type
      conditionalPanel(condition = "input.plotType == 'singleTT'",
        selectInput("tumortype", "Tumor type",
                  c("ACC" = "acc",
                    "BLCA" = "blca",
                    "BRCA" = "brca",
                    "CESC" = "cesc",
                    "CHOL" = "chol",
                    "COAD" = "coad",
                    "DLBC" = "dlbc",
                    "ESCA" = "esca",
                    "GBM" = "gbm",
                    "HNSC" = "hnsc",
                    "KICH" = "kich",
                    "KIRC" = "kirc",
                    "KIRP" = "kirp",
                    "LAML" = "laml",
                    "LGG" = "lgg",
                    "LIHC" = "lihc",
                    "LUAD" = "luad",
                    "LUSC" = "lusc",
                    "MESO" = "meso",
                    "PAAD" = "paad",
                    "PCPG" = "pcpg",
                    "PRAD" = "prad",
                    "READ" = "read",
                    "SARC" = "sarc",
                    "SKCM" = "skcm",
                    "STAD" = "stad",
                    "TGCT" = "tgct",
                    "THCA" = "thca",
                    "THYM" = "thym",
                    "UCEC" = "ucec",
                    "UCS" = "ucs",
                    "UVM" = "uvm",
                    "PanCancer" = "pancancer"))),
      
      conditionalPanel(condition = "input.plotType == 'singleTT'",
        div(checkboxInput(inputId = "advancedSettings",
                          label = "Advanced️︎"),
            align = "right")),
      
      
      # CNA data type
      conditionalPanel(
        condition = "input.plotType == 'singleTT'",
        h3("Copy number"),
        selectInput("cnaType", "Data type",
                  c("GISTIC" = "GISTIC",
                    "Copy number segments" = "CN_Segments"))),
      
      
      # CNA model type
      conditionalPanel(condition = "input.plotType == 'singleTT'",
        selectInput("cnaModel", "Model type",
                  c("Linear model" = "linModel",
                  "Generalized additive model" = "gamMod"))),
      
      # maxK CNA
      conditionalPanel(
        condition = "input.plotType == 'singleTT' && input.cnaModel == 'gamMod' && input.advancedSettings == 1",
        sliderInput("maxK_cna", "Maximum basis functions",
                    min = 3,
                    max = 10,
                    value = 4)),
      
      # Number of PCs to use
      conditionalPanel(condition = "input.plotType == 'singleTT'",
        h3("Methylation"),
        sliderInput("methPCs", "Methylation signatures to use",
                    min = 1,
                    max = 10,
                    value = 4)),
      
      # Meth model type
      conditionalPanel(condition = "input.plotType == 'singleTT'",
        selectInput("methModel", "Model type",
                  c("Generalized additive model" = "gamMod",
                    "Linear model" = "linModel"))),
      
      # maxK Methylation
      conditionalPanel(condition = "input.plotType == 'singleTT' && input.methModel == 'gamMod' && input.advancedSettings == 1",
        sliderInput("maxK_meth", "Maximum basis functions",
                    min = 3,
                    max = 10,
                    value = 4)),
      
      # Meth model type
      conditionalPanel(condition = "input.plotType == 'singleTT'",
        h4("Custom window size"),
        splitLayout(textInput("inputChr", "Chrom"),
                    textInput("inputStart", "Start"),
                    textInput("inputEnd", "End"))),

      actionButton("Submit", label = "Submit")

    ),
    
    ################################################## 
    ## Main panel for displaying outputs
    
    mainPanel(
      
      tabsetPanel(
        
        ## Plotting panel
        tabPanel("Plot",
                 
                 # Text for gene location
                 div(htmlOutput(outputId = "output_geneLocationText")  %>% withSpinner(color = adjustcolor("#4C856E",
                                                                                                           alpha.f = 0.0),
                                                                                       proxy.height = "0px"),
                     align = "center",
                     style = "width: 96%;"),
                 
                 # Plot for chromosome ideogram
                 div(plotOutput(outputId = "ideogram", height = "45px")  %>% withSpinner(color = adjustcolor("#4C856E",
                                                                                                             alpha.f = 0.0),
                                                                                         proxy.height = "0px"),
                 
                 # Break
                 br(),
                 
                 # Top plot (copy number)
                 div(div(plotOutput(outputId = "distPlot", height = "425px")  %>% withSpinner(color="#0dc5c1"),
                     style = "max-width: 800px; display: block; margin-left: auto; margin-right: auto;"),
                  
                 # Bottom plot (methylation)          
                 div(plotOutput(outputId = "distPlotMeth", height = "425px")  %>% withSpinner(color="#0dc5c1"),
                     style = "max-width: 800px; display: block; margin-left: auto; margin-right: auto;"),
                 
                 # Combined plot (in pan-cancer view, if checked off)
                 conditionalPanel(condition = "input.plotType == 'pancancer' && input.showCombined == 1",
                                  div(plotOutput(outputId = "combinedPlot", height = "425px")  %>% withSpinner(color="#0dc5c1"),
                     style = "max-width: 800px; display: block; margin-left: auto; margin-right: auto;")),
                 
                 # Output for R^2 text
                 div(htmlOutput(outputId = "outputRsq_meth")  %>% withSpinner(color = adjustcolor("#4C856E",
                                                                                                  alpha.f = 0.0),
                                                                              proxy.height = "0px"),
                     align = "center"),
                  
                  # Box shadow for all of the above
                 style = "box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19); padding: 10px; max-width: 850px; width: 96%; display: block; margin-left: auto; margin-right: auto;"),
                  
                 # Margins and alignment for the whole div
                 style = "max-width: 800px; width: 96%; display: block; margin-left: auto; margin-right: auto;"),
                 
                 # Break
                 br(),
                 
                 # Download figures
                 div(downloadButton("downloadFig", "Download"),
                     align = "right"),
                 
                 # Break
                 br()),
        
        ## Data panel
        tabPanel("Data",
                 # Table if in pan-cancer mode
                 div(tableOutput("dataTable")  %>% withSpinner(color="#0dc5c1"),
                     align = "left"),
                 
                 # Break
                 br(),
                 
                 # Download button
                 div(downloadButton("downloadData", "Download"),
                     align = "left"),
                 
                 # Break
                 br()),
        
        ## About panel
        tabPanel("About",
                 h2("About"),
                 p(HTML(paste0("PANORAMA (Pan-cancer atlas of transcriptional dependence on DNA methylation and copy number aberrations) is a tool for exploring the associations between gene expression, DNA methylation and copy number in cancer. PANORAMA has two main modes: Pan-cancer, and single model. In the pan-cancer view, it is possible to explore the varying degrees of Expression-Copy number (E-C) and Expression-Methylation (E-M) association across tumor types. A combined model of gene expression as a function of both copy number and methylation (E-CM) is also available. Data in the pan-cancer view are comparable across tumor types. In the single model view, it is possible to perform detailed analyses of E-M and E-C associations in individual tumor types (or to analyze tissue-agnostic/pan-cancer associations). Analyses in the single model view are heavily customizable, and raw data can be downloaded. Model statistics for the single model view are, however, not comparable across tumor types. PANORAMA uses data from ", a(href = "https://www.cell.com/cell/fulltext/S0092-8674(18)30302-7", "The Cancer Genome Atlas"), ". For full methodological details, please read the manuscript associated with PANORAMA, available at bioRxiv."))),
                 br(),
                 p("Please cite any analyses using PANORAMA as follows:"),
                 p("Fougner, C., Höglander, E.K., Lien, T.G., Sørlie, T.S., Nord, S. & Lingjærde, O.C. A pan-cancer atlas of transcriptional dependence on DNA methylation and copy number aberrations.", tags$i("bioRxiv"), "(2020)."),
                 br(),
                 h2("Frequently asked questions"),
                 h4("Are the statistics from models in the pan-cancer view comparable between genes and tumor types?"),
                 p("In the pan-cancer view, model statistics from the same gene should be considered comparable across tumor types. Across genes, model statistics related to copy number should also be considered comparable. Model statistics related to methylation are to an extent comparable across genes, but certain biases introduced by the methylation arrays are unavoidable (e.g. some genes may only be associated with a few CpG probes, while other genes may be associated with several hundred CpG probes). These biases are partially adressed by our PCA-based approach, but should nonetheless be kept in mind when interpreting results."),
                 br(),
                 h4("Are the statistics from models in the pan-cancer view comparable between methylation and CNA?"),
                 p("Methylation and copy number-based model statistics are not directly comparable as different methods are used to model them. Expression-methylation associations are modeled using five covariates (MethSigs) and generalized additive models, whereas expression-methylation associations are modeled using one covariate and linear regression."),
                 br(),
                 h4("Why aren’t all 33 tumor types from The Cancer Genome Atlas shown in the pan-cancer view?"),
                 p("The tumor types in The Cancer Genome Atlas have variable sample numbers. In order to make model statistics comparable across tumor types (in pan-cancer analyses), we randomly downsampled to one hundred tumors, modeled transcriptional associations, and repeated sampling/modeling one hundred times. The median values of these hundred runs were used as the best estimate of model statistics. As a result, all tumor types with fewer than one hundred samples, with all required data, were excluded from pan-cancer analyses. All tumor types (except for ovarian serous cystadenocarcinoma) can be analyzed in the single model view, however statistics here are not comparable across tumor types. Data from Illumina HumanMethylation450 arrays were only available for ten ovarian serous cystadenocarcinoma samples, and the tumor type was therefore fully excluded."),
                 br(),
                 h4("Why do results differ between the pan-cancer view and the single model view?"),
                 p("Results differ between the the pan-cancer view and the single model view due to the downsampling method described above. In the single model view, all available samples are used for modeling."),
                 br(),
                 h4("What units are used in the single model view?"),
                 p(HTML(paste0("Gene expression (RNA-Seq, FPKM-UQ normalized RSEM data) is in the form of ", tags$i("log2(ReadCount + 1)"), ". Copy number (based on Affymetrix SNP6.0 arrays, segments derived using ", a(href = "https://academic.oup.com/biostatistics/article/5/4/557/275197", "circular binary segmentation"), ", made gene centric using Ziggurat Deconstruction in ", a(href = "https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-4-r41", "GISTIC2.0"), ") is in the form of ", tags$i("log2(copynumber/2)"), ", meaning 0 represents a copy number of 2. Methylation data (from HumanMethylation450 arrays) are principal components (MethSigs) derived from all CpGs associated with a gene (by default defined as being within 50 000 bases from coding regions). See the associated article for further details regarding methylation data. MethSigs are directionless and can be mirrored."))),
                 br(),
                 h4("Why can’t I find data for my gene of interest?"),
                 p("Only protein coding genes with data available for all three data levels (gene expression, methylation and copy number) are included in PANORAMA. In the pan-cancer view, only genes associated with at least five CpGs were included, but these genes can be analyzed in the single model view. Genes with variation in expression below a certain level were also excluded."),
                 br(),
                 h4("What do the tumor type abbreviations stand for?"),
                 p("Abbreviations for tumor types can be found", a(href = "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations", "here"), "."),
                 br(),
                 h4("What does the shaded area around the predicted line in single model view represent?"),
                 p("The shaded line represents the 95% confidence interval."),
                 br(),
                 h4("How can I find out which CpGs are used in a methylation analysis?"),
                 p("Raw data from individual probes are not provided if using the standard window size (50 000 bases), as MethSigs are pre-computed for improved performance. Raw probe-level data are however provided when using custom window sizes. To find out which CpGs are included in an analysis using the standard window size, simply enter the gene location (provided above the plot) -/+ 50 000 bases added to the start and end positions."),
                 br(),
                 h4("Where can the source code for PANORAMA be found?"),
                 p(HTML(paste0("The source code for the web application can be found ", a(href = "https://github.com/clfougner/PANORAMA_app", "here"), ". The source code for the analyses underlying results shown in the pan-cancer view can be found ", a(href = "https://github.com/clfougner/PANORAMA", "here"), "."))),
                br(),
                h4("How can the full PANORAMA dataset be downloaded?"),
                p(HTML(paste0("The full PANORAMA dataset can be downloaded from ", a(href = "https://github.com/clfougner/PANORAMA/blob/master/Output/AllModels.zip", "here"), ". Adjusted ", tags$i("R"), tags$sup("2"), " values used for expression-copy number associations are found in the column", tags$i("lm_logE_logC_r2_adj"), ". Adjusted ", tags$i("R"), tags$sup("2"), " values used for expression-methylation associations are found in the column ", tags$i("gam_logE_M_r2_adj"), ". Adjusted ", tags$i("R"), tags$sup("2"), " values used for the combined model are found in the column ", tags$i("methGAM_cnaLM_logE_logC_M_r2_adj"), "."))),
                br(),
                h2("Attributions"),
                p(HTML(paste0("The data underlying all analyses are from ", a(href = "https://www.cell.com/cell/fulltext/S0092-8674(18)30302-7", "The Cancer Genome Atlas pan-cancer dataset"), "."))),
                p(HTML(paste0("The code used for chromosome ideograms is based on code from the ", a(href = "http://www.bioconductor.org/packages/release/bioc/html/copynumber.html", "copynumber R package"), " by ", a(href = "https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-591", "Nilsen et al"), "."))),
                br())
                
      )
    )
  )
)


##################################################
## Server
server <- function(input, output) {
  
   # Read the pan-cancer data
   dat <- readRDS("data/AllDat.rds")
   
   # Read cytoband data
   cytoData <- read.table(file = "data/CytoBands.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
   
   # Read the list of available files
   allFiles <- readRDS("data/AllFiles.rds")
   
   # Read the probe index for custom methylation values
   probeIdx <- readRDS("data/ProbeIndex10kFiles.rds")
   
   # Read the gene position data
   genePositions <- readRDS("data/GenePositions.rds")

   ##################################################
   ## Validate gene in pan-cancer data: with text
   validateCNA_panCan <- reactive({
     
     # Get lower case gene name
     gn <- toupper(input$gene)
     
     # Do nothing if one or zero characters entered
     if(nchar(gn) < 2){
       validate(
         need(nchar(gn) > 1, "")
       )
     }
     
     # Validate
     validate(
       need(gn %in% dat$gene, "\n \n Data not available for the requested gene.")
     )
     
     # Check if gene position available
     validate(
       need(gn %in% genePositions$GeneName, "\n \n Data not available for the requested gene.")
     )
   })
   
   ##################################################
   ## Validate gene in pan-cancer data: no text
   validateMeth_panCan <- reactive({
     
     # Get lower case gene name
     gn <- toupper(input$gene)
     
     # Do nothing if one or zero characters entered
     if(nchar(gn) < 2){
       validate(
         need(nchar(gn) > 1, "")
       )
     }
     
     # Validate
     validate(
       need(gn %in% dat$gene, "")
     )
     
     # Check if gene position available
     validate(
       need(gn %in% genePositions$GeneName, "")
     )
   })
   
   
   ##################################################
   ## Subset the pan-cancer data for a gene
   panCanDat <- reactive({
     
     # If in pan-cancer view
     if(input$plotType == "pancancer"){
       
     # Subset the data
     panCanDat1 <- dat[dat$gene == toupper(input$gene), ]
     
     # Re-order columns
     panCanDat1 <- panCanDat1[, c(5, 2, 3, 4)]
     
     # Return data
     return(panCanDat1)
     }
   })
   
   
   ##################################################
   ## Format the pan-cancer data for copy number
   panCanDatCNA <- reactive({
     
     # If in pan-cancer view
     if(input$plotType == "pancancer"){
       
       # Sort by copy number R^2_adj
       if(input$sortOrder == "individually" | input$sortOrder == "bycna"){
         panCanDatCNA1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_CNA"]), ]    
       }
       
       # Sort by methylation R^2_adj
       if(input$sortOrder == "bymeth"){
         panCanDatCNA1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_Methylation"]), ]    
       }
       
       # Sort by combined R^2_adj
       if(input$sortOrder == "bycombined"){
         panCanDatCNA1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_Combined"]), ]    
       }
       
       # Sort alphabetical
       if(input$sortOrder == "alphabetical"){
         panCanDatCNA1 <- panCanDat()[order(panCanDat()[, "TumorType"], decreasing = TRUE), ]    
       }
     
     # Make a vector of R^2_adj values for copy number
     cd <- panCanDatCNA1[, "Rsq_adj_CNA"]
     
     # Set names to tumor type
     names(cd) <- panCanDatCNA1[, "TumorType"]
     
     # Return R^2_adj values for copy number
     return(cd)
     }
   })
   
   
   ##################################################
   ## Format the pan-can data for methylation
   panCanDatMeth <- reactive({
     
     # If in pan-cancer view
     if(input$plotType == "pancancer"){
       
       # Sort by methylation R^2_adj
       if(input$sortOrder == "individually" | input$sortOrder == "bymeth"){
          panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_Methylation"]), ]
       }
       
       # Sort by copy number R^2_adj
       if(input$sortOrder == "bycna"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_CNA"]), ]
       }
       
       # Sort by combined R^2_adj
       if(input$sortOrder == "bycombined"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_Combined"]), ]
       }
       
       # Sort alphabetical
       if(input$sortOrder == "alphabetical"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "TumorType"], decreasing = TRUE), ]
       }

     # Make a vector of R^2_adj values for methylation
     cd <- panCanDatMeth1[, "Rsq_adj_Methylation"]
     
     # Set names to tumor type
     names(cd) <- panCanDatMeth1[, "TumorType"]

     # Return R^2_adj values for methylation
     return(cd)
     }
   })
   
   
   ##################################################
   ## Format the pan-can data for the combined model
   panCanDatCombined <- reactive({
     
     # If in pan-cancer view
     if(input$plotType == "pancancer"){
       
       # Sort by combined R^2_adj
       if(input$sortOrder == "individually" | input$sortOrder == "bycombined"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_Combined"]), ]
       }
       
       # Sort by copy number R^2_adj
       if(input$sortOrder == "bycna"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_CNA"]), ]
       }
       
       # Sort by methylation R^2_adj
       if(input$sortOrder == "bymeth"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "Rsq_adj_Methylation"]), ]
       }
       
       # Sort alphabetical
       if(input$sortOrder == "alphabetical"){
         panCanDatMeth1 <- panCanDat()[order(panCanDat()[, "TumorType"], decreasing = TRUE), ]
       }
       
       # Make a vector of R^2_adj values for the combined model
       cd <- panCanDatMeth1[, "Rsq_adj_Combined"]
       
       # Set names to tumor type
       names(cd) <- panCanDatMeth1[, "TumorType"]
       
       # Return R^2_adj values for the combined model
       return(cd)
     }
   })
   
   
   ##################################################
   ## Validate single model data: with text
   validateCNA_individual <- reactive({
     
     # Get upper case gene name
     gn <- toupper(input$gene)
     
     # Do nothing if one or zero characters entered
     if(nchar(gn) < 2){
       validate(
         need(nchar(gn) > 1, "")
       )
     }
     
     # If not in custom window size
     if(input$inputChr == "" | input$inputStart == "" | input$inputEnd == ""){
       
       # Build methylation string
       methString <- paste0(input$tumortype,"/Methylation/Precomputed/PCs/", gn, ".rds")
       
       # Validate if methylation string is in the list of available files
       validate(
         need(methString %in% allFiles, "\n \n Data not available for the requested gene.")
       )
       
     }
     
     # If in custom window size
     if(input$inputChr != "" & input$inputStart != "" & input$inputEnd != ""){
       
       # Check if chromosome is correctly formatted
       validate(
         need(toupper(input$inputChr) %in% c(1:22, "X", "Y"),
              "\n \n Chromosome not correctly formatted. Must be a number between 1 and 22, X or Y")
       )
       
       # Get start and end values, remove spaces and commas
       startVal <- gsub(" ", "", input$inputStart)
       startVal <- gsub(",", "", startVal)
       startVal <- as.numeric(startVal)
       
       endVal <- gsub(" ", "", input$inputEnd)
       endVal <- gsub(",", "", endVal)
       endVal <- as.numeric(endVal)
       
       # Check that start and end positions are numbers (as.numeric() above will give NA if not a number)
       validate(
         need(is.na(startVal) == FALSE & is.na(endVal) == FALSE,
              "\n \n Start and end positions must be integers")
       )
       
       # Check that start and end positions are integers
       validate(
         need(startVal %% 1 == 0 & endVal %% 1 == 0,
              "\n \n Start and end positions must be integers")
       )
       
       # Check that start and end positions are positive
       validate(
         need(startVal > 0 & endVal > 0,
              "\n \n Start and end positions must be positive values")
       )
       
       # Check that start is before end
       validate(
         need(startVal<= endVal,
              "\n \n End position must be greater than or equal to the start position")
       )
       
       # Check that the chromosomal range is within the accepted limits
       validate(
         need((endVal - startVal) <= 10000000,
              "\n \n Maximum permitted window size: 10 000 000 base pairs. Please reduce window size.")
       )
       
       # Check if there are any CpGs in the requested position
       currCpGs <- probeIdx[probeIdx$Chrom == toupper(input$inputChr), ]
       currCpGs <- currCpGs[currCpGs$Pos >= startVal & currCpGs$Pos <= endVal, ]
       validate(
         need(nrow(currCpGs) > 0,
              "\n \n No CpGs available within the window.")
       )
       
     }
     
     # Build expression string
     exprString <- paste0(input$tumortype, "/Expression/", gn, ".rds")
     
     # Build CNA string
     cnaString <- paste0(input$tumortype, "/CNA/", input$cnaType, "/", gn, ".rds")
     
     # Validate expression file
     validate(
       need(exprString %in% allFiles, "\n \n Data not available for the requested gene.")
     )
     
     # Validate CNA file
     validate(
       need(cnaString %in% allFiles, "\n \n Data not available for the requested gene.")
     )
     
     # Check if gene position available
     validate(
       need(gn %in% genePositions$GeneName, "\n \n Data not available for the requested gene.")
     )
     
   })
   
   
   ##################################################
   ## Validate single model data: no text
   validateOther_individual <- reactive({
     
     # Get upper case gene name
     gn <- toupper(input$gene)
     
     # Do nothing if one or zero characters entered
     if(nchar(gn) < 2){
       validate(
         need(nchar(gn) > 1, "")
       )
     }
     
     # If not in custom window size
     if(input$inputChr == "" | input$inputStart == "" | input$inputEnd == ""){
       
       # Build methylation string
       methString <- paste0(input$tumortype,"/Methylation/Precomputed/PCs/", gn, ".rds")
       
       # Validate if methylation string is in the list of available files
       validate(
         need(methString %in% allFiles, "")
       )
     }
     
     # If in custom window size
     if(input$inputChr != "" & input$inputStart != "" & input$inputEnd != ""){
       
       # Check if chromosome is correctly formatted
       validate(
         need(toupper(input$inputChr) %in% c(1:22, "X", "Y"),
              "")
       )
       
       # Get start and end values
       startVal <- gsub(" ", "", input$inputStart)
       startVal <- gsub(",", "", startVal)
       startVal <- as.numeric(startVal)
       
       endVal <- gsub(" ", "", input$inputEnd)
       endVal <- gsub(",", "", endVal)
       endVal <- as.numeric(endVal)
       
       # Check that start and end positions are numbers (as.numeric() above will give NA if not a number)
       validate(
         need(is.na(startVal) == FALSE & is.na(endVal) == FALSE,
              "")
       )
       
       # Check that start and end positions are integers
       validate(
         need(startVal %% 1 == 0 & endVal %% 1 == 0,
              "")
       )
       
       # Check that start and end positions are positive
       validate(
         need(startVal > 0 & endVal > 0,
              "")
       )
       
       # Check that start is before end
       validate(
         need(startVal<= endVal,
              "")
       )
       
       # Check that the chromosomal range is within the accepted limits
       validate(
         need((endVal - startVal) <= 10000000,
              "")
       )
       
       # Check if there are any CpGs in the requested position
       currCpGs <- probeIdx[probeIdx$Chrom == toupper(input$inputChr), ]
       currCpGs <- currCpGs[currCpGs$Pos >= startVal & currCpGs$Pos <= endVal, ]
       validate(
         need(nrow(currCpGs) > 0,
              "")
       )
     }
     
     # Build expression string
     exprString <- paste0(input$tumortype, "/Expression/", gn, ".rds")
     
     # Build CNA string
     cnaString <- paste0(input$tumortype, "/CNA/", input$cnaType, "/", gn, ".rds")
     
     # Validate expression file
     validate(
       need(exprString %in% allFiles, "")
     )
     
     # Validate CNA file
     validate(
       need(cnaString %in% allFiles, "")
     )
     
     # Check if gene position available
     validate(
       need(gn %in% genePositions$GeneName, "")
     )
   })
   
   ##################################################
   ## Validate MethPCs and K against sample size for tumor type
   illegalCombinations <- reactive({
     
     # Tumor type
     tt <- input$tumortype
     
     # Make string for GAM using selected number of PCs and max basis functions (k)
     currentCombination <- paste0(input$methPCs, ", ", input$maxK_meth)
     
     # Create "val" object == FALSE, i.e. default to no error for validation
     val <- FALSE
     
     # Set "illegal" to empty if no disallowed combinations (i.e. there are enough samples for all possible model parameters)
     illegal <- c()
     
     ## Check parameters based on sample numbers in tumor type
     
     if(tt == "acc"){
       # MethPCs, K. Number of tumors = 76
       illegal <- c("10, 10",
                    "10, 9",
                    "9, 10")
     }
     
     if(tt == "chol"){
       # MethPCs, K, number of tumors = 36
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "10, 6",
                    "10, 5",
                    "10, 4",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "9, 7",
                    "9, 6",
                    "9, 5",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "8, 7",
                    "8, 6",
                    "7, 10",
                    "7, 9",
                    "7, 8",
                    "7, 7",
                    "6, 10",
                    "6, 9",
                    "6, 8",
                    "6, 7",
                    "5, 10",
                    "5, 9",
                    "4, 10")
     }
     
     if(tt == "dlbc"){
       # MethPCs, K, number of tumors = 47
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "10, 6",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "9, 7",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "8, 7",
                    "7, 10",
                    "7, 9",
                    "7, 8",
                    "6, 10",
                    "6, 9")
     }
     
     if(tt == "gbm"){
       # MethPCs, K, number of tumors = 54
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "9, 7",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "7, 10",
                    "7, 9",
                    "6, 10")
     }
     
     if(tt == "kich"){
       # MethPCs, K, number of tumors = 61
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "8, 10",
                    "8, 9",
                    "7, 10")
     }
     
     if(tt == "meso"){
       # MethPCs, K, number of tumors = 83
       illegal <- c("10, 10")
     }
     
     if(tt == "read"){
       # MethPCs, K, number of tumors = 87
       illegal <- c("10, 10")
     }
     
     if(tt == "ucs"){
       # MethPCs, K, number of tumors = 56
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "7, 10",
                    "7, 9")
     }
     
     if(tt == "uvm"){
       # MethPCs, K, number of tumors = 80
       illegal <- c("10, 10",
                    "10, 9",
                    "9, 10")
     }
     
     # Is the current combination allowed?
     val <- currentCombination %in% illegal
     
     # Validate
     validate(
       need(val != TRUE, "Too many MethSigs and/or basis functions relative to sample number for the tumor type. Please reduce MethSigs and/or basis functions.")
     )
   })
   
   # No text
   illegalCombinationsNoText <- reactive({
     
     # Tumor type
     tt <- input$tumortype
     
     # Make string for the current combination of PCs and gams
     currentCombination <- paste0(input$methPCs, ", ", input$maxK_meth)
     
     # Set "illegal" to empty if no disallowed combinations (i.e. there are enough samples for all possible model parameters)
     illegal <- c()
     
     # Create "val" object == FALSE, i.e. default to no error for validation
     val <- FALSE
     
     ## Check parameters based on sample numbers in tumor type
     
     if(tt == "acc"){
       # MethPCs, K. Number of tumors = 76
       illegal <- c("10, 10",
                    "10, 9",
                    "9, 10")
     }
     
     if(tt == "chol"){
       # MethPCs, K, number of tumors = 36
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "10, 6",
                    "10, 5",
                    "10, 4",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "9, 7",
                    "9, 6",
                    "9, 5",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "8, 7",
                    "8, 6",
                    "7, 10",
                    "7, 9",
                    "7, 8",
                    "7, 7",
                    "6, 10",
                    "6, 9",
                    "6, 8",
                    "6, 7",
                    "5, 10",
                    "5, 9",
                    "4, 10")
     }
     
     if(tt == "dlbc"){
       # MethPCs, K, number of tumors = 47
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "10, 6",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "9, 7",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "8, 7",
                    "7, 10",
                    "7, 9",
                    "7, 8",
                    "6, 10",
                    "6, 9")
     }
     
     if(tt == "gbm"){
       # MethPCs, K, number of tumors = 54
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "9, 7",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "7, 10",
                    "7, 9",
                    "6, 10")
     }
     
     if(tt == "kich"){
       # MethPCs, K, number of tumors = 61
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "8, 10",
                    "8, 9",
                    "7, 10")
     }
     
     if(tt == "meso"){
       # MethPCs, K, number of tumors = 83
       illegal <- c("10, 10")
     }
     
     if(tt == "read"){
       # MethPCs, K, number of tumors = 87
       illegal <- c("10, 10")
     }
     
     if(tt == "ucs"){
       # MethPCs, K, number of tumors = 56
       illegal <- c("10, 10",
                    "10, 9",
                    "10, 8",
                    "10, 7",
                    "9, 10",
                    "9, 9",
                    "9, 8",
                    "8, 10",
                    "8, 9",
                    "8, 8",
                    "7, 10",
                    "7, 9")
     }
     
     if(tt == "uvm"){
       # MethPCs, K, number of tumors = 80
       illegal <- c("10, 10",
                    "10, 9",
                    "9, 10")
     }
     
     # Is the current combination allowed?
     val <- currentCombination %in% illegal
     
     # Validate
     validate(
       need(val != TRUE, "")
     )
   })
   
   ##################################################
   ## Get methylation data for individual model view
   getMethDat <- reactive({
  
  # Get upper case gene name
  gn <- toupper(input$gene)
     
   # If not in custom window size
   if(input$inputChr == "" | input$inputStart == "" | input$inputEnd == ""){
     
     # Read MethSigs/PCs
     methDat <- readRDS(paste0("/data/panemec/PanEMEC/",
                               input$tumortype,
                               "/Methylation/Precomputed/PCs/",
                               gn,
                               ".rds"))
     
     # Rename to MethSigs
     colnames(methDat) <- paste0("MethSig", 1:ncol(methDat))
     
     # Get variance captured
     vars <- readRDS(paste0("/data/panemec/PanEMEC/",
                            input$tumortype,
                            "/Methylation/Precomputed/VarianceCaptured/",
                            gn,
                            ".rds"))
     
     # Round to 3 digits
     vars <- round(vars, digits = 3)
     
     # Return methylation data
     return(list(methDat, vars))
     
   }
   
   # If in custom window size
   if(input$inputChr != "" & input$inputStart != "" & input$inputEnd != ""){

     # Get start and end values
     startVal <- gsub(" ", "", input$inputStart)
     startVal <- gsub(",", "", startVal)
     startVal <- as.numeric(startVal)
     
     endVal <- gsub(" ", "", input$inputEnd)
     endVal <- gsub(",", "", endVal)
     endVal <- as.numeric(endVal)
     
     # Which files are needed?
     currCpGs <- probeIdx[probeIdx$Chrom == toupper(input$inputChr), ]
     currCpGs <- currCpGs[currCpGs$Pos >= startVal & currCpGs$Pos <= endVal, ]
     
     files <- unique(currCpGs$FileIndex)
     
     # Read the first file
     dat <- readRDS(paste0("/data/panemec/PanEMEC/",
                           input$tumortype,
                           "/Methylation/BetaValues/",
                           files[1],
                           ".rds"))
     
     # If more than one file is needed
     if(length(files) > 1){
       # For each file
       for(i in files){
         # Read file
         currentDat <- readRDS(paste0("/data/panemec/PanEMEC/",
                                      input$tumortype,
                                      "/Methylation/BetaValues/",
                                      i,
                                      ".rds"))
         # CBind to the dat object
         dat <- cbind(dat, currentDat)
       }
     }
     
     # Subset the required CpGs
     dat <- dat[, currCpGs$Probe]

     # Run PCA
     currentPCA <- prcomp(x = dat, scale = FALSE, center = TRUE)
     
     # Save to methDat object
     methDat <- cbind(currentPCA$x, dat)
     
     # Rename to MethSigs
     colnames(methDat)[grep(x = colnames(methDat), pattern = "PC")] <- paste0("MethSig", 1:length(grep(x = colnames(methDat), pattern = "PC")))
     
     # Get variance
     vars <- summary(currentPCA)[["importance"]]["Proportion of Variance", ]
     
     # Round to three digits
     vars <- round(vars, digits = 3)
     
     # Return methylation data
     return(list(methDat, vars))
     
   }
})

   
   ##################################################
   ## Get data for individual models: gene expression, copy number and methylation. Methylation calls the above function
   individualData <- reactive({
     if(input$plotType == "singleTT"){
       
       # Get upper case gene name
       gn <- toupper(input$gene)
       
       # Read gene expression data
       exprDat <- readRDS(paste0("/data/panemec/PanEMEC/",
                                 input$tumortype,
                                 "/Expression/",
                                 gn,
                                 ".rds"))
       
      # Get methylation data
       methDat <- getMethDat()

       
       # Read copy number data
       cnDat <- readRDS(paste0("/data/panemec/PanEMEC/",
                               input$tumortype,
                               "/CNA/",
                                input$cnaType,
                                "/",
                                gn,
                                ".rds"))

       # Bind all of the above to a single data frame
       dataForGam <- cbind(exprDat, cnDat, methDat[[1]])
       
       # Remove NAs
       dataForGam <- dataForGam[complete.cases(dataForGam), ]
       
       ## Add jitter to CNA data
       # Max and min value for CNA
       cna_cap <- 3.657
       
       if(input$cnaType == "GISTIC"){
         cna_cap_low <- -1.293
       }
       
       if(input$cnaType == "CN_Segments"){
         cna_cap_low <- -14.2877
       }
       
       
       # Add jitter to maximum value
       if(max(dataForGam[, "CopyNumber"], na.rm = TRUE) == cna_cap){
          dataForGam[dataForGam[, "CopyNumber"] == cna_cap, "CopyNumber"] <- jitter(dataForGam[dataForGam[, "CopyNumber"] == cna_cap, "CopyNumber"])
       }
       
       # Add jitter to minimum value
       if(min(dataForGam[, "CopyNumber"], na.rm = TRUE) == cna_cap_low){
          dataForGam[dataForGam[, "CopyNumber"] == cna_cap_low, "CopyNumber"] <- jitter(dataForGam[dataForGam[, "CopyNumber"] == cna_cap_low, "CopyNumber"])
       }
       
       # Return data for individual model
       return(data.frame(dataForGam))
     }
   })
   

   ##################################################
   ## Make pan-cancer copy number plot
   makePanCan_CNA_plot <- function(){
     
     # Validate data
     validateCNA_panCan()
     
     # Make plot
     p <- barplot(panCanDatCNA(),
             xlab = parse(text = '~italic(R[Adjusted]^"2")'),
             horiz = TRUE,
             main = "Expression - Copy number",
             xlim = c(0, 1),
             las = 1,
             font.main = 1,
             col = "#046C9A")
     
     # Return plot
     return(p)
   }
   
   
   ##################################################
   ## Make single model copy number plot
   makeSingle_CNA_plot <- function(){
     
     # Check if data is available
     validateCNA_individual()
     
     # Use non-linear or linear model?
     # If linear model
     if(input$cnaModel == "linModel"){
       
       # Fit linear model
       gm <- lm(Expression ~ CopyNumber,
                data = individualData())
       
       # Get rsq_adj
       rsq <- round(summary(gm)$r.squared, digits = 3)
       
       # Predict from linear regression
       sef <- predict(gm, se.fit = TRUE, interval = "confidence")
       
       # Bind data together
       d <- cbind(pred = sef$fit[, "fit"],
                  se_up = sef$fit[, "lwr"],
                  se_dn = sef$fit[, "upr"],
                  cn = individualData()[, "CopyNumber"])
     }
     
     # If GAM
     if(input$cnaModel == "gamMod"){
       
       # Run model
       gm <- gam(formula = Expression ~ s(CopyNumber, k = input$maxK_cna),
                 data = individualData(),
                 method = "REML")
       
       # Get rsq_adj
       rsq <- round(summary(gm)$dev.expl, digits = 3)
       
       # Predict from generalized additive model
       sef <- predict.gam(gm, se.fit = TRUE)
       
       # Bind data together
       d <- cbind(pred = sef$fit,
                  se_up = sef$fit + sef$se.fit * 1.96,
                  se_dn = sef$fit - sef$se.fit * 1.96,
                  cn = individualData()[, "CopyNumber"])
     }
     
     # Order data
     d <- d[order(d[, "cn"]), ]
     
     # Plot E-C association
     plot(Expression ~ CopyNumber,
          data = individualData(),
          main = bquote(italic(R^2)==.(rsq)),
          font.main = 1,
          pch = 16,
          cex = 0.5,
          col = "#046C9A")
     
     # Add confidence interval
     polygon(x = c(d[, "cn"],
                   d[order(d[, "cn"], decreasing = TRUE), "cn"]),
             y = c(d[, "se_up"],
                   d[order(d[, "cn"], decreasing = TRUE), "se_dn"]),
             col = adjustcolor("#046C9A", alpha.f = 0.2),
             border = NA)
     
     # Add predicted line (white background line for contrast against points)
     lines(pred ~ cn,
           data = d,
           lwd = 4,
           col = "white")
     
     # Add predicted line
     lines(pred ~ cn,
           data = d,
           lwd = 2,
           col = "#046C9A")
    
   }
   
   
   ##################################################
   ## Make PanCancer methylation plot
   makePanCan_Meth_plot <- function(){
     
     # Validate data
     validateMeth_panCan()
     
     # Make plot
     barplot(panCanDatMeth(),
             horiz = TRUE,
             xlab = parse(text = '~italic(R[Adjusted]^"2")'),
             main = "Expression - Methylation",
             xlim = c(0, 1),
             las = 1,
             font.main = 1,
             col = "#4C856E")
     }
   
   
   ##################################################
   ## Make single model methylation plot
   makeSingle_Meth_plot <- function(){
     
     # Validate data
     validateOther_individual()
     
     # Validate k and MethPCs
     illegalCombinations()
     
     # Get columns (MethSigs) available for methylation data
     pcsAvailable <- length(grep(x = colnames(individualData()), pattern = "MethSig"))
     
     # If there are fewer MethSigs available than selected (input$MethPCs), use those available
     if(pcsAvailable < input$methPCs){
       methPCs <- pcsAvailable
     }
     
     # If there are more MethSigs available than selected (input$MethPCs), use the selected number of MethSigs
     if(pcsAvailable >= input$methPCs){
       methPCs <- input$methPCs
     }
     
     # Number of rows and columns for plots, depending on number of MethSigs used
     if(methPCs == 1){r <- 1; c <- 1}
     if(methPCs == 2){r <- 1; c <- 2}
     if(methPCs == 3){r <- 1; c <- 3}
     if(methPCs == 4){r <- 2; c <- 2}
     if(methPCs == 5 | methPCs == 6){r <- 2; c <- 3}
     if(methPCs == 7 | methPCs == 8){r <- 2; c <- 4}
     if(methPCs == 9 | methPCs == 10){r <- 2; c <- 5}
     
     # Set plotting parameters (rows and columns)
     par(mfrow = c(r, c))
     
     # For each MethSig
     for(i in 1:methPCs){
       
       # Use non-linear or linear model?
       if(input$methModel == "linModel"){
         
         ## If linear model
         # Generate formula
         frm <- reformulate(termlabels = paste0("MethSig", i),
                            response = "Expression")
         
         # Run model
         gm <- lm(formula = frm,
                  data = individualData())
         
         # Get rsq_adj
         rsq <- round(summary(gm)$r.squared, digits = 3)
         
         
         # Predict
         sef <- predict(gm, se.fit = TRUE, interval = "confidence")
         
         # Bind data together
         d <- cbind(pred = sef$fit[, "fit"],
                    se_up = sef$fit[, "lwr"],
                    se_dn = sef$fit[, "upr"],
                    individualData()[, paste0("MethSig", i)])
         
         # Rename to MethSig
         colnames(d)[4] <- paste0("MethSig", i)
       }
       
       ## If generalized additive model
       if(input$methModel == "gamMod"){
         
         # Generate formula
         frm <- reformulate(termlabels = paste0("s(MethSig",
                                                i,
                                                ", k = ",
                                                input$maxK_meth, ")"),
                            response = "Expression")
         
         # Run model
         gm <- gam(formula = frm,
                   data = individualData(),
                   method = "REML")
         
         # Get rsq_adj
         rsq <- round(summary(gm)$dev.expl, digits = 3)
         
         
         # Predict
         sef <- predict.gam(gm, se.fit = TRUE)
         
         # Bind data together
         d <- cbind(pred = sef$fit,
                    se_up = sef$fit + sef$se.fit * 1.96,
                    se_dn = sef$fit - sef$se.fit * 1.96,
                    individualData()[, paste0("MethSig", i)])
         
         # Rename to MethSig
         colnames(d)[4] <- paste0("MethSig", i)
       }
       
       # Order data
       d <- d[order(d[, paste0("MethSig", i)]), ]
       
       # Plot E-M model
       plot(reformulate(termlabels = paste0("MethSig", i),
                        response = "Expression"),
            data = individualData(),
            pch = 16,
            cex = 0.5,
            main = bquote(atop(Proportion~of~variance ==  .(getMethDat()[[2]][i]), italic(R^2)==.(rsq))),
            font.main = 1,
            col = "#4C856E")
       
       # Add confidence interval
       polygon(x = c(d[, paste0("MethSig", i)],
                     d[order(d[, paste0("MethSig", i)], decreasing = TRUE), paste0("MethSig", i)]),
               y = c(d[, "se_up"],
                     d[order(d[, paste0("MethSig", i)], decreasing = TRUE), "se_dn"]),
               col = adjustcolor("#4C856E", alpha.f = 0.2),
               border = NA)
       
       # Add predicted line (white background line for contrast against points)
       lines(reformulate(termlabels = paste0("MethSig", i),
                         response = "pred"),
             data = d,
             lwd = 4,
             col = "white")
       
       # Add predicted line
       lines(reformulate(termlabels = paste0("MethSig", i),
                         response = "pred"),
             data = d,
             lwd = 2,
             col = "#4C856E")
     }
   }
   
   
   ##################################################
   ## Make PanCancer combined plot
   makePanCan_Combined_plot <- function(){
     
     # Validate data
     validateMeth_panCan()
     
     # Make plot
     barplot(panCanDatCombined(),
             horiz = TRUE,
             xlab = parse(text = '~italic(R[Adjusted]^"2")'),
             main = "Expression - Combined",
             xlim = c(0, 1),
             las = 1,
             font.main = 1,
             col = "#C93312")
   }
   
   
   ##################################################
   ## Get gene position data for top text
   geneLocationDat <- reactive({
     
     # If in single model view
     if(input$plotType == "singleTT"){
       
       # Validate data
       validateOther_individual()
       
       # Validate k and MethPCs
       illegalCombinationsNoText()
     }
     
     # If in pan-cancer view
     if(input$plotType == "pancancer"){
       
       # Validate data
       validateMeth_panCan()
       
     }
     
     # Get upper case gene name
     gn <- toupper(input$gene)
     
     # Get gene position data
     positionDat <- genePositions[genePositions$GeneName == gn, ]
     
     # Get gene chromosome
     geneChrom <- as.character(positionDat$Chrom[1])
     
     # Get gene start position
     geneStart <- positionDat$Start[1]
     
     # Format start position (spaces every third digit)
     geneStartFormatted <- formatC(as.numeric(positionDat$Start[1]),
                                   big.mark = " ",
                                   big.interval = 3, 
                                   format = "d",
                                   flag = "0",
                                   width = nchar(as.numeric(positionDat$Start[1])))
     
     # Get gene end position
     geneEnd <- positionDat$End[1]
     
     # Format end position (spaces every third digit)
     geneEndFormatted <- formatC(as.numeric(positionDat$End[1]),
                                 big.mark = " ",
                                 big.interval = 3, 
                                 format = "d",
                                 flag = "0",
                                 width = nchar(as.numeric(positionDat$End[1])))
     
     # Get karyotype band
     geneBand <- as.character(positionDat$KaryotypeBand)
     
     # Add chromosome to karyotype band
     karyoFormatted <- paste0(geneChrom, geneBand)
     
     # Make named vector of all relevant gene position data
     result <- c("GeneName" = gn,
                 "Chrom" = geneChrom,
                 "Start" = geneStart,
                 "StartFormatted" = geneStartFormatted,
                 "End" = geneEnd,
                 "EndFormatted" = geneEndFormatted,
                 "GeneBand" = geneBand,
                 "KaryoFormatted" = karyoFormatted)
     
     # Return the above vector
     return(result)
   })
   
   
   ##################################################
   ## Make chromosome ideogram
   ## This code is adapted from copynumber package by Nilsen et al.
   ## http://www.bioconductor.org/packages/release/bioc/html/copynumber.html
   ## https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-591
   
   plotIdeogram <- function(chrom, cex = 0.6, cyto.data, cyto.unit = "bp", unit, posStart, posEnd){
     
      # Subset cytoband data
     chrom.cytoband <- cyto.data[cyto.data[, 1] == paste("chr", chrom, sep = ""), ]

     # Start positions
     cyto.start <- chrom.cytoband[, 2]
     
     # End position
     cyto.end <- chrom.cytoband[, 3]
     
     # Left side x-coordinate
     xleft <- cyto.start
     
     # Right side x-coordinate
     xright <- cyto.end
     
     # Number of cytobands
     n <- length(xleft)
     
     # Length of the current chromosome
     chrom.length <- xright[n] - xleft[1]
     
     # Stain data
     stain <- chrom.cytoband[, 5]
     
     # List of stain types
     sep.stain <- c("gpos", "gneg", "acen", "gvar", "stalk")
     
     # Find rows with each stain type
     g <- sapply(sep.stain, grep, x = stain, fixed = TRUE)
     
     # Centromere
     centromere <- g$acen
     
     # Stalk
     stalk <- g$stalk

     # Set colors
     col <- rep("", n)
     co <- "#262F46"
     col[stain == "gneg"] <- "#fdfdfd"
     col[stain == "gpos100"] <- adjustcolor(co, alpha.f = 1)
     col[stain == "gpos75"] <- adjustcolor(co, alpha.f = 0.75)
     col[stain == "gpos50"] <- adjustcolor(co, alpha.f = 0.5)
     col[stain == "gpos25"] <- adjustcolor(co, alpha.f = 0.25)
     col[stain == "stalk"] <- adjustcolor(co, alpha.f = 0.1)
     col[stain == "gvar"] <- adjustcolor(co, alpha.f = 0.75)
     col[stain == "acen"] <- "#FDD262"
     
     # Shading/lines
     density <- rep(NA, n)
     
     # Line angles
     angle <- rep(45, n)
     
     # Line density
     density[stain == "gvar"] <- 15
     
     # Lower y-limit
     ylow <- 0
     
     # High y-limit
     yhigh <- 1
     
     # Make plot
     plot(x = c(0, max(xright)),
          y = c(ylow, yhigh), 
          type = "n",
          axes = FALSE,
          xlab = "",
          ylab = "",
          xlim = c(0, max(xright)),
          ylim = c(0, 1),
          xaxs = "i")
     
     # Rectangles
     skip.rect <- c(1, centromere, n, stalk)
     
     # Make rectangles
     rect(xleft[-skip.rect],
          rep(ylow, n - length(skip.rect)),
          xright[-skip.rect],
          rep(yhigh, n - length(skip.rect)),
          col = col[-skip.rect],
          border = "black",
          density = density[-skip.rect],
          angle = angle[-skip.rect])
     
     # Round edges at ideogram start, stop and at centromere
     draw.roundEdge(start = xleft[1],
                    stop = xright[1],
                    y0 = ylow,
                    y1 = yhigh,
                    col = col[1],
                    bow = "left",
                    density = density[1],
                    angle = angle[1],
                    chrom.length = chrom.length)
     
     draw.roundEdge(start = xleft[centromere[1]],
                    stop = xright[centromere[1]],
                    y0 = ylow,
                    y1 = yhigh,
                    col = col[centromere[1]],
                    bow = "right",
                    density = density[centromere[1]],
                    angle = angle[centromere[1]],
                    lwd = 1,
                    chrom.length = chrom.length)
     
     draw.roundEdge(start = xleft[centromere[2]],
                    stop = xright[centromere[2]],
                    y0 = ylow,
                    y1 = yhigh,
                    col = col[centromere[2]],
                    bow = "left",
                    density = density[centromere[2]],
                    angle = angle[centromere[2]],
                    lwd = 1,
                    chrom.length = chrom.length)
     
     draw.roundEdge(start = xleft[n],
                    stop = xright[n],
                    y0 = ylow,
                    y1 = yhigh,
                    col = col[n],
                    bow = "right",
                    density = density[n],
                    angle = angle[n],
                    chrom.length = chrom.length)
     
     # Draw stalk-segment
     if(length(stalk) > 0){
       for(i in 1:length(stalk)){
         drawStalk(xleft[stalk[i]], xright[stalk[i]], ylow, yhigh, col = col[stalk[i]])
       }
     }
     
     # Add dot for chromosome position
     stripchart(x ~ y,
                data = data.frame(x = as.numeric(posStart), y = 0),
                at = 0.5,
                bg = "#FD6467",
                col = "black",
                pch = 23,
                cex = 2,
                add = TRUE)
     
   }
   
   ## Draw round edge function for chromosome ideogram
   # This code is adapted from copynumber package by Nilsen et al.
   # http://www.bioconductor.org/packages/release/bioc/html/copynumber.html
   # https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-591
   draw.roundEdge <- function(start, stop, y0, y1, col, bow, density = NA, angle = 45, lwd = 1, chrom.length){
     
     # Y points in round edge
     f <- rep(0, 0)
     f[1] <- 0.001
     i = 1
     half <- y0 + (y1 - y0)/2
     
     while(f[i] < half){
       f[i+1] <- f[i] * 1.3
       i <- i+1
     }
     
     f <- f[-length(f)]
     
     Y <- c(y1, y1, y1 - f, half, y0 + rev(f), y0, y0)
     
     # X points in roundedge
     cyto.length <- stop - start
     
     share <- cyto.length/chrom.length
     if(share > 0.2){
       
       # To create bow in end of chromosome 24
       share <- 0.2
     }
     
     if(bow == "left"){
       
       round.start <- start + cyto.length * (1 - share)^20
       
       x <- seq(round.start, start, length.out = (length(f) + 2))
       revx <- rev(x[-length(x)])
       x <- c(x, revx)
       X <- c(stop, x, stop)
     } else {
       if(bow == "right"){
         round.start <- stop - cyto.length * (1 - share)^20
         x <- seq(round.start, stop, length.out = (length(f) + 2))
         revx <- rev(x[-length(x)])
         x <- c(x, revx)
         X <- c(start, x, start)
       }
     }
     
     polygon(x = X,
             y = Y,
             col = col,
             border = "black",
             density = density, 
             angle = angle,
             lwd = lwd)
   }
   
   ## Draw stalk for chromosome ideogram
   # This code is adapted from copynumber package by Nilsen et al.
   # http://www.bioconductor.org/packages/release/bioc/html/copynumber.html
   # https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-591
   drawStalk <- function(start, stop, y0, y1, col){
     
     size <- stop - start
     x1 <- c(start, start + size/3, stop - size/3, stop)
     x2 <- rev(x1)
     x <- c(x1, x2)
     y_1 <- c(y0, y0 + 0.25, y0 + 0.25, y0)
     y_2 <- c(y1 , y1 - 0.25, y1 - 0.25, y1)
     y <- c(y_1, y_2)
     polygon(x = x, y = y, col = col)
     
   }
   
   
   ##################################################
   ## Render chromosome ideogram
   output$ideogram <- renderPlot({
      
      # Get gene location data
      dat <- geneLocationDat()
      
      # Set margins
      par(mar = c(0, 1, 1, 1))
      
      # Return plot
      return(plotIdeogram(chrom = dat["Chrom"], cyto.data = cytoData, posStart = dat["Start"], posEnd = dat["End"]))
   })
  
  ##################################################
  ## Render copy number plot
  output$distPlot <- renderPlot({
    
    # If in pan-cancer view
    if(input$plotType == "pancancer"){
        
      # Set margins
      par(mar = c(4, 5, 5, 1)) # Bottom, left, top, right
      
      # Return pan-cancer plot
      return(makePanCan_CNA_plot())
    }
    
    # If in single model view
    if(input$plotType == "singleTT"){
      
      # Return single model plot
      return(makeSingle_CNA_plot())
    }
  })
  
  ##################################################
  ## Render methylation plot
  output$distPlotMeth <- renderPlot({

    # If in pan-cancer view
    if(input$plotType == "pancancer"){
      
      # Set margins
      par(mar = c(4, 5, 5, 1)) # Bottom, left, top, right
      
      # Return pan-cancer plot
      return(makePanCan_Meth_plot())
     
    }
    
    # If in single model view
    if(input$plotType == "singleTT"){
      
      # Return single model plot
      return(makeSingle_Meth_plot())
    }
  })
   
   ##################################################
   ## Render combined plot
   output$combinedPlot <- renderPlot({
     
     # If in pan-cancer view and checked off on "Show combined"
     if(input$plotType == "pancancer" & input$showCombined == 1){
       
       # Set margins
       par(mar = c(4, 5, 5, 1))
       
       # Return pan-cancer plot
       return(makePanCan_Combined_plot())
     }
   })
   
  
  ##################################################
  ## Render table on "Data" panel (only if in pan-cancer view)
  output$dataTable <- renderTable({
    
    # If in single tumor type view
    if(input$plotType == "singleTT"){
      
      # Validate data so that error message is returned, even though a table isn't rendered
      validateCNA_individual()
    }
    
    # If in pan-cancer view
    if(input$plotType == "pancancer"){
      
      # Validata data
      validateCNA_panCan()
      
      # Return data for table
      return(panCanDat())
    }
  })
  
   
  ##################################################
  ## Download figure
  output$downloadFig <- downloadHandler(
    
    # Set file name
    filename = function() {
      
      # If in single model view
      if(input$plotType == "singleTT"){
        
        # Make file name string
        fn <- paste0(toupper(input$tumortype), "_", toupper(input$gene), "_singleModels.pdf")
      }
      
      # If in pan-cancer view
      if(input$plotType == "pancancer"){
        
        # Make file name string
        fn <- paste0(toupper(input$gene), "_panCanModels.pdf")
      }
      
      # Return file name string
      return(fn)
    },
    
    # Set content of download
    content = function(file) {
      
      # If in single model view
      if(input$plotType == "singleTT"){
        
        # Get columns available for methylation data
        pcsAvailable <- ncol(individualData()) - 2
        
        # If fewer MethSigs are available than the number selected, use the available MethSigs
        if(pcsAvailable < input$methPCs){
          methPCs <- pcsAvailable
        }
        
        # If enough MethSigs are available, use the number of MethSigs selected
        if(pcsAvailable >= input$methPCs){
          methPCs <- input$methPCs
        }
        
        # Set figure width and height based on the number of MethSigs used
        if(methPCs == 1){w <- 4; h = 4}
        if(methPCs == 2){w <- 5; h = 4}
        if(methPCs == 3){w <- 6; h = 4}
        if(methPCs == 4){w <- 5; h = 4}
        if(methPCs %in% c(5, 6)){w <- 6; h = 4}
        if(methPCs %in% c(7, 8)){w <- 7; h = 4}
        if(methPCs %in% c(9, 10)){w <- 8; h = 4}
        
        # Make PDF
        pdf(file, height = h, width = w, pointsize = 8, useDingbats = FALSE)
          
          # Copy number plot
          makeSingle_CNA_plot()
          
          # Methylation plot
          makeSingle_Meth_plot()
          
        # End PDF
        dev.off()
      }
      
      # If in pan-cancer view
      if(input$plotType == "pancancer"){
        
        # If combined model not shown
        if(input$showCombined == 0){
          
          # Make PDF
          pdf(file, height = 4, width = 5, pointsize = 8, useDingbats = FALSE)
            
            # Set number of columns in figure
            par(mfrow = c(1, 2))
            
            # Copy number plot
            makePanCan_CNA_plot()
            
            # Methylation plot
            makePanCan_Meth_plot()
          
          # End PDF
          dev.off()  
        }
  
        # If combined model shown
        if(input$showCombined == 1){
          
          # Make pDF
          pdf(file, height = 4, width = 7, useDingbats = FALSE)
          
          # Set number of columns in figure
            par(mfrow = c(1, 3))
            
            # Make copy number plot
            makePanCan_CNA_plot()
            
            # Make methylation plot
            makePanCan_Meth_plot()
            
            # Make combined plot
            makePanCan_Combined_plot()
            
          # End plot
          dev.off()  
        }
      }
    }
  )
  
  
  ##################################################
  ## Download data
  output$downloadData <- downloadHandler(
    
    # Set file name
    filename = function() {
      
      # If in single model view
      if(input$plotType == "singleTT"){
        
        # Make file name string
        fn <- paste0(toupper(input$tumortype), "_", toupper(input$gene), "_data.txt")
      }
      
      # If in pan-cancer view
      if(input$plotType == "pancancer"){
      
        # Make file name string
        fn <- paste0(toupper(input$gene), "_models.txt")
      }
      
      # Return file name
      return(fn)
    },
    
    # Set content of download
    content = function(file) {
      
      # If in single model view
      if(input$plotType == "singleTT"){
        
        # Validate data
        validateCNA_individual()

        # Get data
        dat <- individualData()
        
        # Make a column for sample IDs
        dat$SampleID <- rownames(dat)
        
        # Re-order columns
        dat <- dat[, c(ncol(dat), 1:(ncol(dat) - 1))]
        
        # Write table
        write.table(dat, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
      
      # If in pan-cancer view
      if(input$plotType == "pancancer"){
        
        # Validata data
        validateCNA_panCan()
        
        # Write table
        write.table(panCanDat(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  )
  
  ##################################################
  ## Render the bottom text
  output$outputRsq_meth <- renderUI({
    
    # If in single model view
    if(input$plotType == "singleTT"){
      
      # Validate data
      validateOther_individual()
      
      # Validate k and MethPCs
      illegalCombinationsNoText()
      
      # Get columns available for methylation data
      pcsAvailable <- length(grep(x = colnames(individualData()), pattern = "MethSig"))
      
      # If fewer MethSigs are available than selected (input$methPCs), use the number of MethSigs available
      if(pcsAvailable < input$methPCs){
        methPCs <- pcsAvailable
      }
      
      # If enough MethSigs are available, use the number of MethSigs selected
      if(pcsAvailable >= input$methPCs){
        methPCs <- input$methPCs
      }
      
      ## Rsquared for the full methylation model with all MethSigs
      # If linear model
      if(input$methModel == "linModel"){
        
        # Generate formula string
        frm <- reformulate(termlabels = paste0("MethSig", 1:methPCs),
                           response = "Expression")
        
        # Run model
        gm <- lm(formula = frm,
                 data = individualData())
        
        # Get rsq_adj
        rsq_m <- round(summary(gm)$r.squared, digits = 3)
        
        ## Combined model
        # If copy number == linear regression
        if(input$cnaModel == "linModel"){
          
          # Generate formula string
          frm_cm <- reformulate(termlabels = c(paste0("MethSig",
                                                      1:methPCs),
                                               "CopyNumber"),
                                response = "Expression")
          
          # Run model
          gm_cm <- lm(formula = frm_cm,
                      data = individualData())
          
          # Get Rsq_adj
          rsq_cm <- round(summary(gm_cm)$r.squared, digits = 3)
          
        }
        
        # If copy number == generalized additive model
        if(input$cnaModel == "gamMod"){
          
          # Generate formula string
          frm_cm <- reformulate(termlabels = c(paste0("MethSig",
                                                      1:methPCs),
                                                paste0("s(CopyNumber, k =",
                                                      input$maxK_cna,
                                                      ")")),
                                response = "Expression")
          
          # Run model
          gm_cm <- gam(formula = frm_cm,
                      data = individualData(),
                      method = "REML")
          
          # Get Rsq_adj
          rsq_cm <- round(summary(gm_cm)$dev.expl, digits = 3)
        }  
      }
      
      # If generalized additive model
      if(input$methModel == "gamMod"){
        
        # Generate formula string
        frm <- reformulate(termlabels = paste0("s(MethSig",
                                               1:methPCs,
                                               ", k = ",
                                               input$maxK_meth,
                                               ")"),
                           response = "Expression")
        
        # Run model
        gm <- gam(formula = frm,
                  data = individualData(),
                  method = "REML")
        
        # Get rsq_adj
        rsq_m <- round(summary(gm)$dev.expl, digits = 3)
        
        ## Combined model
        # If copy number == linear regression
        if(input$cnaModel == "linModel"){
          
          # Generate formula string
          frm_cm <- reformulate(termlabels = c(paste0("s(MethSig",
                                                      1:methPCs,
                                                      ", k = ",
                                                      input$maxK_meth,
                                                      ")"),
                                               "CopyNumber"),
                                response = "Expression")
          
          # Run model
          gm_cm <- gam(formula = frm_cm,
                   data = individualData(),
                   method = "REML")
          
          # Get Rsq_adj
          rsq_cm <- round(summary(gm_cm)$dev.expl, digits = 3)
          
        }
        
        # If copy number == generalized additive model
        if(input$cnaModel == "gamMod"){
          
          # Generate formula string
          frm_cm <- reformulate(termlabels = c(paste0("s(MethSig",
                                                      1:methPCs,
                                                      ", k = ",
                                                      input$maxK_meth,
                                                      ")"),
                                               paste0("s(CopyNumber, k =",
                                                      input$maxK_cna,
                                                      ")")),
                                response = "Expression")
          
          # Run model
          gm_cm <- gam(formula = frm_cm,
                       data = individualData(),
                       method = "REML")
          
          # Get Rsq_adj
          rsq_cm <- round(summary(gm_cm)$dev.expl, digits = 3)
        }
      }
      
      # More than or equal to one principal component?
      if(methPCs > 1){pluralPCs <- "s"}
      if(methPCs == 1){pluralPCs <- ""}
      
      # Get gene position and karyotype band
      gn <- toupper(input$gene)
      positionDat <- genePositions[genePositions$GeneName == gn, ]
      geneChrom <- positionDat$Chrom[1]
      geneStart <- positionDat$Start[1]
      geneEnd <- positionDat$End[1]
      geneBand <- positionDat$KaryotypeBand
      geneSentence <- ""
      
      # Get the cumulative variance captured by MethSigs
      cumVarCaptured <- sum(getMethDat()[[2]][1:methPCs])
      
      # If cumVarCaptured is greater than 1 (due to rounding), set to 1
      if(cumVarCaptured > 1){cumVarCaptured <- 1}
      
      # Generate bottom string
      txt <- paste0("<br> <b>Adjusted <i>R<sup>2</sup></i> for the full methylation model : ",
                   rsq_m,
                   "<br>Adjusted <i>R<sup>2</sup></i> for the combined methylation/copy number model : ",
                   rsq_cm,
                  "</b> <br>Cumulative proportion of variance captured by ",
                  methPCs,
                  " methylation signature",
                  pluralPCs,
                  " : ",
                  cumVarCaptured,
                  "<br> Number of samples: ",
                  nrow(individualData()))
      
      # If fewer MethSigs are available than selected
      if(pcsAvailable < input$methPCs){
        
        # Make string informing that fewer MethSigs are available than selected
        txt <- paste0(txt,
                      "<br> More methylation signatures selected than are available. Using ",
                      pcsAvailable,
                      " methylation signature",
                      pluralPCs, ".")
      }
      
      # Render text
      HTML(txt)
    }
  })
   

  ##################################################
  ## Render the gene position: text
  output$output_geneLocationText <- renderUI({
      
      # Get gene location data
      dat <- geneLocationDat()
      
      # Make gene location string
      locationSentence <- paste0("<br><i>",
                                 dat["GeneName"],
                                 "</i>",
                                 " located on ",
                                 dat["KaryoFormatted"],
                                 " (position: ",
                                 dat["StartFormatted"],
                                 " - ",
                                 dat["EndFormatted"],
                                 ")")

      # Render text
      HTML(locationSentence)
    
  })
}


shinyApp(ui = ui, server = server)
