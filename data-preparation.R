# Libraries required
library(dplyr)
library(HGNChelper)
library(tidyr)
library(here)
library(stringr)

# Prepare the data

ingest_data_default <- function(param) {
  
  data <- read.csv(paste0(param$dataset, ".csv"))
  subset <- data %>% dplyr::select(
    GeneID,
    starts_with("logFC"),
    starts_with("adjpval"),
    starts_with("dataset")
  )
  subset$dataset <- param$dataset
  colnames(subset) <- c("GeneID", "logFC", "padj","dataset")
  return(subset)
}

load_data <- function() {
  database <- read.csv(file.path('metadata.csv'))
  print(typeof(database))
  datasets <- list()

  for(entry in 1:nrow(database)) {
    tmp <- ingest_data_default(database[entry, ])
    datasets[[database[entry, "dataset"]]] <- tmp
    #datasets <- append(datasets, var)
   # print(head(datasets))
  }
str(datasets)
return(datasets)
}

load_data()
