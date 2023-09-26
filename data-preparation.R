# Libraries required
library(dplyr)
library(HGNChelper)
library(tidyr)
library(here)
library(stringr)

# Prepare the data

update_genenames <- function(){
#  date <- as.character(Sys.Date())
#  filename <- paste0("HGNC_genenames_", date, ".RData")
  filename <- paste0("HGNC_genenames_newest", ".RData")
  genenames_newest <- getCurrentHumanMap()
  save(genenames_newest, file = filename) 
}

update_genenames()

ingest_data_default <- function(param) {
  data <- read.csv(paste0(param$dataset, ".csv"))
  subset <- data %>% dplyr::select(
    GeneID,
    starts_with("logFC"),
    starts_with("adjpval"),
    starts_with("dataset")
  )
  subset$GeneID <- sapply(strsplit(subset$GeneID,";"), `[`, 1) # in case there are multiple IDs, take the first one
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
  
  load("HGNC_genenames_newest.RData")
  for (entry in 1:nrow(database)) {
    genenames_update <- checkGeneSymbols(database[[entry]]["geneID"], species = "human", map = genenames_newest, unmapped.as.na = FALSE)
    database[[entry]]["geneID"] <- genenames_update$Suggested.Symbol #replace the original with updated IDs
  }
return(datasets)
}

load_data()
