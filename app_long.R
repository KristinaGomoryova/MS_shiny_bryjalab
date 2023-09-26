library(here)
library(shiny)
library(shinythemes)
library(shinymanager)
library(plotly)
library(ggplot2)
library(DT)
library(dplyr)
library(gprofiler2)
library(rentrez)
library(mygene)

###################################    Credentials   #######################################
credentials <- read.delim(here('credentials.txt'))

#####################################   Datasets   #########################################
# Get the datasets
# Stored in the /data folder
data_2006.RNF43 <- read.csv(here('data', 'shinyapp_2006-RNF43.csv'))
data_2434.2660 <- read.csv(here('data', 'shinyapp_2434-2660.csv'))
data_2006.RNF122 <- read.csv(here('data', 'shinyapp_2006-RNF122.csv'))
data_2006.RNFT2 <- read.csv(here('data', 'shinyapp_2006-RNFT2.csv'))
data_2999 <- read.csv(here('data', 'shinyapp_2999.csv'))
data_3019 <- read.csv(here('data', 'shinyapp_3019.csv'))
data_4310.CK1 <- read.csv(here('data', 'shinyapp_4310-CK1s.csv'))
data_4310.FZD <- read.csv(here('data', 'shinyapp_4310-FZDs.csv'))
data_4453.DDA <- read.csv(here('data', 'shinyapp_4453-DDA.csv'))
data_4453.DIA <- read.csv(here('data', 'shinyapp_4453-DIA.csv'))

datasets.long <- read.delim(here("datasets_long.txt"))
select.protein <- sort(unique(datasets.long$GeneID))

################################ Data description tab  #####################################
# Data description tab
dataDescription <- data.frame(
  dataset = c("2006_RNF43",
              "2006_RNF122",
              "2006_RNFT2",
              "2434_2660",
              "2999",
              "3019",
              "3676",
              "4310_CK1s",
              "4310_FZDs",
              "4453_DDA",
              "4453_DIA"),
  acquisition = c("DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DDA",
                  "DIA"),
  responsible_person = c("Tomek",
                         "Tomek",
                         "Tomek",
                         "Tomek",
                         "Marek",
                         "Tomek/PetaP",
                         "Vendy",
                         "Tomas",
                         "Tomas",
                         "Tomas",
                         "Tomas"),
  PRIDE = c("PXD020478",
            "-",
            "-",
            "PXD039259",
            "PXD034237, PXD033548",
            "-",
            "PXD041751",
            "-",
            "-",
            "-",
            "-"),
  location = c("Orders\\2006_RNF43\\final_table.xlsx",
               "Orders\\2006_RNF122\\RNF122_for_BK\\updated_KNIME_tables-BirA-removed.xlsx",
               "Orders\\2006_RNFT2\\2006_RNFT2_publication\\results-table_20221109.xlsx",
               "Orders\\2434_2660\\final-for-publication\\DEP_processing_BH.xlsx",
               "Orders\\2999\\2999_KNIME_20230802.txt",
               "Orders\\3019\\version_13\\analysis\\outputs\\01_DEP-processing\\01_DEP-results-table.csv",
               "Orders\\3676\\project_1_UC-vs-qEV\\version_06\\2023_Vyhlidalova_JEV_submission\\outputs\\03_data-processed\\03_data-processed.csv",
               "Orders\\4310\\CK1s_no-wnt_batch-corrected_FINAL.html",
               "Orders\\4310\\Fzd_batch-corrected_FINAL.html",
               "Orders\\4453\\4453_interaction-limma_results.xlsx",
               "Orders\\3906_4453_CK1s\\4453\\4453_report.xlsx"),
  analysis = c("KNIME",
               "KNIME",
               "KNIME",
               "DEP",
               "KNIME",
               "DEP",
               "R",
               "DEP",
               "DEP",
               "DEP",
               "KNIME"),
  note = c("",
           "qualitative hits reported separately, check the source file on server",
           "",
           "",
           "log2, norm; check source for imp",
           "",
           "no stats done",
           "batch correction",
           "batch correction",
           "interaction limma",
           ""),
  download_MM = c(
    '<a href="https://kristinagomoryova.shinyapps.io/ms_shiny_bryjalab/MM_2006_RNF43.docx"> "MM 2006 RNF43"</a>',
    '<a href="https://kristinagomoryova.shinyapps.io/ms_shiny_bryjalab/MM_2434_2660.docx"> "MM 2434-2660"</a>',
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    ""))

description_page <- tabPanel(
  title = "Description",
  mainPanel(
    DT::dataTableOutput("description")      
  )
)

################################  Explore dataset tab  #####################################

exploreDataset_page <- tabPanel(
  title = "Explore dataset",
  titlePanel = "Explore dataset",
  fluidRow(
    column(3,
           selectInput("dataset", 
                       "Choose a dataset:",
                       choices = c("data_2006.RNF43",
                                   "data_2006.RNF122",
                                   "data_2006.RNFT2",
                                   "data_2434.2660",
                                   "data_2999",
                                   "data_3019",
                                   "data_3676",
                                   "data_4310.CK1",
                                   "data_4310.FZD",
                                   "data_4453.DDA",
                                   "data_4453.DIA",
                                   "test"),
                       selected = "data_2006.RNF43"),
           
           selectInput("logFC", "Log2 fold change (logFC):", c()),
           selectInput("pvalue", "Adjusted p-value (adjpval):", c()),
           selectInput("annotation", "Annotation:", c()),
           selectInput("geneOntology",
                       "Run gene ontology?",
                       choices = c("none", "on upregulated", "on downregulated"),
                       selected = "none")),
    
    # Show a plot of the generated distribution
    column(9,   
           h2("Statistical analysis of the dataset"),
           DT::dataTableOutput("table")
    )),
  fluidRow(
    column(9, offset = 3,
           h2("Volcano plot"),
           plotlyOutput("volcanoplot")),
    column(9, offset = 3,       
           h2("Selected points from volcano plot"),
           DT::dataTableOutput("selected_data")
    ),
    column(9, offset = 3,
           h2("Gene ontology: gProfiler2"),
           plotlyOutput("gprofilerPlot")))
)

###################################  Find protein tab  #####################################

findProtein_page <- tabPanel(
  title = "Find protein",
  titlePanel = "Find protein",
  fluidRow(
    column(2,
           selectInput("findProtein", 
                       "Find protein:",
                       choices = select.protein,
                       selected = "DVL3",
                       multiple = TRUE)),
    column(10,
           h2("Protein was found in the following datasets: "),
           DT::dataTableOutput("retrievedData")      
    )),
  fluidRow(
    column(2),
    column(10, 
           h2("Entrez database information"),
           textOutput("entrezDescription"))
  )
)

###################################  How to use this app tab  #####################################

howTo_page <- tabPanel(
  title = "How to use this app",
  titlePanel = "How to use this app",
  hr(),
  h2("Welcome to the database of mass spectrometry experiments in Bryjalab!"),
  br(),
  p("Here you will find more details on how to run this app and what you can do in individual tabs: "),
  br(),
  tags$ul(
    tags$li(tags$b("Explore dataset"), " - in this tab you can explore the datasets we produced, plot the volcano plot and perform Gene ontology analysis using gprofiler2"),
    tags$li(tags$b("Find protein"), " - in this tab you can find protein of your interest across datasets and whether it was up/downregulated in any contrast"),
    tags$li(tags$b("Description"), " - in this tab you can find description of particular proteomic experiment")
  ),
  br(),
  p("What to pay attention to: "),
  tags$ul(
    tags$li("To display a volcano plot, all three arguments: logFC, adjusted pvalue and annotation MUST be selected"),
    tags$li("To inspect selected proteins from the volcano plot, use the 'Box Select' tool "),
    tags$li("To run gene ontology, Annotation MUST be set to gene names or majority protein IDs")
  ),
  p("Gene names were updated on 2023-09-20.")
)

#########################################  UI  #############################################
ui <- navbarPage(
  title = "MS datasets in Bryjalab",
  theme = shinytheme('united'),
  verbatimTextOutput("auth_output"),
  exploreDataset_page,
  findProtein_page,
  description_page,
  howTo_page
)

ui <- secure_app(ui)
#######################################  server  ###########################################

server <- function(input, output, session){
  
  options(shiny.maxRequestSize=10*1024^2)
  
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  
  # Description tab output
  output$description <-  DT::renderDataTable(dataDescription, escape = FALSE)
  
  # Explore dataset tab output
  datasetInput <- eventReactive(input$dataset,{
    get(input$dataset)
  })
  
  observeEvent(input$dataset, {
    req(datasetInput())
    choices <- names(datasetInput())
    updateSelectInput(session, "logFC", choices = choices, selected = choices[1])
    updateSelectInput(session, "pvalue", choices = choices, selected = choices[1])
    updateSelectInput(session, "annotation", choices = choices, selected = choices[1])
  }, ignoreNULL = FALSE)
  
  output$table <- DT::renderDT({
    datatable(datasetInput(), filter="top",
              extensions = 'Buttons',
              options = list(
                buttons = c('copy', 'csv', 'excel'),
                dom = 'Bfrtip',
                scrollY = 300,
                scroller = TRUE,
                scrollX=TRUE
              ),
              class = "display"
    )
  })
  
  data.filtered <- reactive({
    req(datasetInput(), input$logFC, input$pvalue, input$annotation)
    if (input$logFC == "X" & input$pvalue == "X" & input$annotation == "X") {
      return()
    }
    selected_df <- datasetInput() %>% dplyr::select(input$logFC, input$pvalue, input$annotation, X)
    colnames(selected_df) <- c("logFC", "pvalue", "annotation", "X")
    
    selected_df <- selected_df %>%
      mutate(col = case_when(
        pvalue < 0.05 & logFC > 1 ~ "up",
        pvalue < 0.05 & logFC < -1 ~ "down",
        TRUE ~ "ns"
      ))
    
    selected_df
  })
  
  
  output$volcanoplot <- renderPlotly({
    req(data.filtered)
    df = data.filtered()
    p <- df %>%
      ggplot(aes(x = logFC, y = -log10(pvalue), color = col, 
                 customdata = X))+
      geom_point(aes(text = annotation))+
      theme_minimal() + 
      scale_color_manual(values=c("#56B4E9", "#999999", "#E69F00")) +
      labs(x = "logFC",
           y = "-log10(adjusted p-value)")
    
    fig <- ggplotly(p, tooltip = c("text"))
    event_register(fig, "plotly_selected")
    fig
    
  })
  
  # Selected points from the volcano plot
  output$selected_data <- DT::renderDataTable({
    plotly_event_data <- event_data(event = "plotly_selected", priority = "event")
    req(plotly_event_data)
    filter(datasetInput(), X %in% plotly_event_data$customdata)
  })
  
  # Gene ontology analysis
  
  output$gprofilerPlot <- renderPlotly({
    
    df = data.filtered()
    req(input$geneOntology)
    req(input$annotation)
    
    if (input$geneOntology == "on upregulated") {
      if (length(df$annotation[df$col == "up"] > 0)){
        gostres <- gost(query = df$annotation[df$col == "up"], organism = "hsapiens")
        gostplot(gostres, capped = FALSE, interactive = TRUE)} else { print("Zero upregulated proteins")}
    } else if (input$geneOntology == "on downregulated") {
      if (length(df$annotation[df$col == "down"] > 0)){
        gostres <- gost(query = df$annotation[df$col == "down"], organism = "hsapiens")
        gostplot(gostres, capped = FALSE, interactive = TRUE)} else { print("Zero downregulated proteins")}
    } else {
    }
  })
  
  
  # Tab for finding proteins
  
  output$retrievedData <- DT::renderDataTable({
    protein_to_found <-  req(input$findProtein)
    req(datasets.long)
    
    datasets.long <- datasets.long %>%
      mutate(change = case_when(
        logFC > 1 & padj < 0.05 ~ "upregulated",
        logFC < -1 & padj < 0.05 ~ "downregulated",
        TRUE ~ "ns"
      ))
    
    datatable(datasets.long[datasets.long$GeneID %in% protein_to_found, ],
              filter="top",
              extensions = 'Buttons',
              options = list(
                buttons = c('copy', 'csv', 'excel'),
                dom = 'Bfrtip',
                scrollY = 300,
                scroller = TRUE,
                scrollX=TRUE
              ),
              class = "display")
  })
  
  proteinEntrezInfo <- reactive({
    protein_to_find <-  req(input$findProtein)
    protein_to_find <- as.character(protein_to_find)
    entrezSymbol <- queryMany(protein_to_find, scopes="symbol", fields="entrezgene", species="human")
    eg <- entrezSymbol$entrezgene
    entrezInfo <- entrez_summary(db = "gene", id=eg)
    return(entrezInfo$summary)
  })
  
  output$entrezDescription <- renderText({
    paste0(as.character(proteinEntrezInfo()))
  })
  
}

# Run the application   
shinyApp(ui = ui, server = server)
