# Libraries required
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

# Insert the functions
source("data-preparation.R")

data_long <- load_data()
data_wide <- load_data_wide()

#select.protein <- sort(unique(datasets.long$GeneID))
select.protein <- sort(Reduce(c, sapply(data_wide,function(x) x[1]) ))

################################ Data description tab  #####################################
# Data description tab
dataDescription <- read.csv("database/metadata.csv")
dataDescription$type <- NULL

description_page <- tabPanel(
  title = "Description",
  mainPanel(
    DT::dataTableOutput("description")      
  )
)

names(data_wide)

################################  Explore dataset tab  #####################################

exploreDataset_page <- tabPanel(
  title = "Explore dataset",
  titlePanel = "Explore dataset",
  fluidRow(
    column(3,
           selectInput("dataset", 
                       "Choose a dataset:",
                       choices = names(data_wide),
                       selected = "2006_RNF43"),
           
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
  exploreDataset_page,
  findProtein_page,
  description_page,
  howTo_page
)

#######################################  server  ###########################################

server <- function(input, output, session){
  
  # Description tab output
  output$description <-  DT::renderDataTable(dataDescription, escape = FALSE)
  
  # Explore dataset tab output
  datasetInput <- eventReactive(input$dataset,{
    #get(input$dataset)
    data_wide[[input$dataset]]
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
   
    if(length(colnames(selected_df)) == 4) {
      colnames(selected_df) <- c("logFC", "pvalue", "annotation", "X")
      selected_df <- selected_df %>%
        mutate(col = case_when(
          pvalue < 0.05 & logFC > 1 ~ "up",
          pvalue < 0.05 & logFC < -1 ~ "down",
          TRUE ~ "ns"
        ))
    } else {
      selected_df <- data.frame(
        logFC=numeric(0),
        pvalue=numeric(0),
        annotation=character(0),
        X=numeric(0),
        col=character(0)
      )
    }

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
    protein_to_find <-  req(input$findProtein)
  #  req(datasets.long)
    req(data_long)
    
    data_long[["3676"]] <- NULL
    data_long_df <- bind_rows(data_long)
    
    data_long_df <- data_long_df %>%
      mutate(change = case_when(
        logFC > 1 & padj < 0.05 ~ "upregulated",
        logFC < -1 & padj < 0.05 ~ "downregulated",
        TRUE ~ "ns"
      ))
    
    datatable(data_long_df[data_long_df$GeneID %in% protein_to_find, ],
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
