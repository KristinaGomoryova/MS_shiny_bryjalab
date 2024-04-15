###############################################################################################################
#################### Update gene names in all datasets to the most current version ############################
###############################################################################################################

# Libraries required
library(HGNChelper)
library(dplyr)
library(here)

# Get the newest genenames update
updateCurrentMap <- function() {
  date <- as.character(Sys.Date())
  filename <- paste0("HGNC_genenames_", date, ".RData")
  genenames_newest <- getCurrentHumanMap()
  save(genenames_newest, file = filename)
}

# Update gene names
updateGeneNames <- function() {
  # Get the input datasets in the database directory:
  tmp.path <- "C:/Users/kika2/Documents/shinyapp_20240411/MS_shiny_bryjalab/database"
#  filenames <- list.files(path=getwd(), pattern = ".csv", full.names = TRUE)
  filenames <- list.files(path=tmp.path, pattern = ".csv", full.names = TRUE)
  
  # Get the newest genename
  updateCurrentMap()
  genenames_list <- file.info(list.files(pattern = "HGNC_", path = ".", full.names = TRUE))
  HGNC <- gsub("/MS_shiny_bryjalab/database/", "", rownames(genenames_list)[which.max(genenames_list$mtime)])
  load(HGNC)
  
  for (i in seq_along(filenames)) {
    if (grepl("metadata.csv", filenames[i])) {
      next
    }

    tmp <- read.csv(filenames[i], row.names=NULL)
    tmp.update <- checkGeneSymbols(
      tmp$GeneID,
      species = "human",
      map = genenames_newest,
      unmapped.as.na = FALSE)
    tmp$GeneID <- tmp.update$Suggested.Symbol
    write.csv(tmp, file=filenames[i], row.names=FALSE)
  }
}

# Run the function
updateGeneNames()

