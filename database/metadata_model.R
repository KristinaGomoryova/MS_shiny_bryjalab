# Libraries required
library(here)

# Metadata table

metadata <- data.frame(
  dataset = c("2006_RNF43"),
  acquisition = c("DDA"),
  responsible_person = c("Tomek"),
  PRIDE = c("PXD020478"),
  analysis = c("KNIME"),
  note = c(""),
  type = c("default"))


# Write the metadata
write.csv(metadata, here('database', 'metadata_model.csv'))
