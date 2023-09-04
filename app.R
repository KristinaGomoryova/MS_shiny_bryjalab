library(here)
library(shiny)
library(shinythemes)
library(shinymanager)
library(plotly)
library(ggplot2)
library(DT)
library(dplyr)
library(gprofiler2)


###################################    Credentials   #######################################
credentials <- read.delim(here('credentials.txt'))

#####################################   Datasets   #########################################
# Get the datasets
# Stored in the /data folder
data_2006.RNF43 <- read.csv(here('data', 'shinyapp_2006-RNF43.csv'))
data_2434.2660 <- read.csv(here('data', 'shinyapp_2434-2660.csv'))

datasets.long <- read.delim(here("datasets_long.txt"))
select.protein <- unique(datasets.long$GeneID)

################################ Data description tab  #####################################
# Data description tab
dataDescription <- data.frame(
  dataset = c("2006_RNF43",
              "2434_2660"),
  acquisition = c("DDA",
                  "DDA"),
  responsible_person = c("Tomek",
                         "Tomek"),
  PRIDE = c("PXD020478",
            "PXD039259"),
  location = c("Orders\\2006_RNF43\\final_table.xlsx",
               "Orders\\2434_2660\\final-for-publication\\DEP_processing_BH.xlsx"),
  analysis = c("KNIME",
               "DEP"),
  note = c("",
           ""),
  download_MM = c(
    '<a href="http://127.0.0.1:5535/MM_2006_RNF43.docx"> "MM 2006 RNF43"</a>',
    '<a href="http://127.0.0.1:5535/MM_2434_2660.docx"> "MM 2434-2660"</a>')
)

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
                                   "data_2434.2660"),
                       selected = "data_2006.RNF43"),
           
           selectInput("logFC", "Log fold change:", c()),
           selectInput("pvalue", "Adjusted p-value:", c()),
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
                     selected = "DVL3")),
    column(10,
         h2("Protein was found in the following datasets: "),
         DT::dataTableOutput("retrievedData")      
  ))
)

#########################################  UI  #############################################
ui <- navbarPage(
  title = "MS datasets in Bryjalab",
  theme = shinytheme('united'),
  verbatimTextOutput("auth_output"),
  exploreDataset_page,
  findProtein_page,
  description_page
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
    selected_df <- datasetInput() %>% select(input$logFC, input$pvalue, input$annotation, X)
    colnames(selected_df) <- c("logFC", "pvalue", "annotation", "X")
    
    selected_df <- selected_df %>%
      mutate(col = case_when(
        pvalue < 0.05 & logFC > 1 ~ "up",
        pvalue < 0.05 & logFC < -1 ~ "down",
        TRUE ~ "ns"
      ))
    print(selected_df)
    selected_df
  })
  

  output$volcanoplot <- renderPlotly({
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

}

# Run the application   
shinyApp(ui = ui, server = server)
