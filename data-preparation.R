# Libraries required
library(dplyr)
library(HGNChelper)
library(tidyr)
library(here)
library(stringr)

# Prepare the data
METADATA_PATH <- here('database', 'metadata.csv')

update_genenames <- function(){
  #  date <- as.character(Sys.Date())
  #  filename <- paste0("HGNC_genenames_", date, ".RData")
  filename <- paste0("HGNC_genenames_newest", ".RData")
  genenames_newest <- getCurrentHumanMap()
  save(genenames_newest, file = filename) 
}

#update_genenames()

ingest_data_default <- function(param) {
  data <- read.csv(here("database", paste0(param$dataset, ".csv")))
  subset <- data %>% dplyr::select(
    GeneID, X,
    starts_with("logFC"),
    starts_with("adjpval")
  )
  
  colnames(subset) <- gsub("adjpval_", "adjpvalX", colnames(subset))
  colnames(subset) <- gsub("logFC_", "logFCX", colnames(subset))
  subset$GeneID <- sapply(strsplit(subset$GeneID,";"), `[`, 1) # in case there are multiple IDs, take the first one
  
  subset <-  pivot_longer(subset, 
                          cols = !c(GeneID,X),
                          names_to = c(".value", "contrast"),
                          names_sep = "X")
  subset <- subset[ , c(1,3,4,2)]
  
  subset$dataset <- param$dataset
  colnames(subset) <- c("GeneID", "logFC", "padj", "contrast", "dataset")
  subset$contrast <- stringr::str_replace_all(subset$contrast, "\\.", "_vs_") #Escape the dot!
  return(subset)
}

ingest_data_noStats <- function(param) {
  data <- read.csv(here("database", paste0(param$dataset, ".csv")))
  data$GeneID <- sapply(strsplit(data$GeneID,";"), `[`, 1) # in case there are multiple IDs, take the first one
  data$dataset <- param$dataset
  data <- data[, c(2, 1, 3:ncol(data))]
  return(data)
}

load_data <- function() {
  database <- read.csv(METADATA_PATH)
  datasets <- list()
  
  for(entry in 1:nrow(database)) {
    if (database[entry, "type"] == "default"){
      datasets[[database[entry, "dataset"]]] <- ingest_data_default(database[entry, ])
    } else if (database[entry, "type"] == "noStats") {
      datasets[[database[entry, "dataset"]]] <- ingest_data_noStats(database[entry, ])
    }
  } 
  return(datasets)
}

#load_data()

ingest_data_wide <- function(param) {
  data <- read.csv(here("database", paste0(param$dataset, ".csv")))
  subset <- data %>% dplyr::select(
    GeneID, X,
    starts_with("logFC"),
    starts_with("adjpval")
  )
  subset$GeneID <- sapply(strsplit(subset$GeneID,";"), `[`, 1) # in case there are multiple IDs, take the first one
  return(subset)
}


load_data_wide <- function() {
  database <- read.csv(METADATA_PATH)
  datasets <- list()
  for(entry in 1:nrow(database)) {
      datasets[[database[entry, "dataset"]]] <- ingest_data_wide(database[entry, ])
    }
  return(datasets)
}

#load_data_wide()
