# Shiny application for sharing mass spectrometry data in Bryjalab

This repository contains instructions on how to run a shiny application, developed to facilitate the inspection of mass spectrometry (MS) datasets produced in [Bryjalab](https://github.com/bryjalab).

The shiny application is composed of 3 modules:

- 'Explore dataset' module allows visual inspection of the data in an interactive table, visualization of selected contrasts in a volcano plot and optional gene ontology analysis on either upregulated or downregulated proteins
- 'Find protein' module allows searching for protein of interest across all available datasets and its associated statistics
- 'Description' module provides an information how the dataset was produced in our lab and on the proteomics core facility

The application was written in R, version 4.3.1, using [Shiny library](https://shiny.posit.co/), version 1.7.5. The web version of application can be accessed only from the MUNI IP address range for the security reasons - for this reason, also data (containing statistical results of contrasts of interest) are stored within a private repository. 

The application can be run also locally, using the model dataset, which is a part of already published manuscript.

# How to run the application

**Running the model data**: 

Currently the repository contains two metadata.R scripts: `metadata.R` containing information about datasets produced in Bryjalab (which are private), and `metadata_model.R`, which is intended to be used to run the model data. In case own data would be added, change the metadata information accordingly.

## Running the shiny app using RStudio

1. Clone this github repository locally
2. Copy the datasets to be displayed into the `/database` folder (be sure you added them also to the `database/metadata.R`). Model data are already in the `/database` folder.
3. Run the `app.R` script

## Running the shiny app using Docker

1. Clone this github repository locally 
2. Copy the datasets to be displayed into the `/database` folder (be sure you added them also to the `database/metadata.R`). Model data are already in the `/database` folder.
3. Mount the data to the docker image, using a similar command as e.g.:
`docker run -p 3838:3838 -v /home/ubuntu/MS_shiny_bryjalab/database/:/MS_shiny_bryjalab/database/ kristinagomoryova/app`
4. Run the shiny application, either as http://localhost:3838, or http://*public_IP*:3838

## Adding a new dataset
New dataset should be added to the `/database` folder. It is preferred, that upon addition of a new dataset, all datasets would get updated gene names in order to synchronize them. A model dataset, `2006_RNF43.csv`, related to publication by [Radaszkiewicz et al.(2021)](https://elifesciences.org/articles/65759) is provided.

Gene names update can be done running the updateGeneNames.R script, which, however, requires installation of R, and specific packages, `dplyr`, `here` and `HGNCHelper`. 

This can be done using commands:
`sudo apt install r-base-core`
`sudo apt-get install r-cran-dplyr`
`R`
`install.packages("HGNCHelper")`

The gene names update script can be run in the terminal using `Rscript updateGeneNames.R`

Don't forget to update the metadata.R and also metadata.csv produced by the R script to include new dataset! 

The new dataset will be displayed automatically upon refreshing the application in the browser.

