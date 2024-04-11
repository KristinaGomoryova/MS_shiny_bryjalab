
# Use an official R runtime as a parent image
FROM rocker/shiny:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
libcurl4-gnutls-dev \
libssl-dev \
libbz2-dev \
liblzma-dev

# install Rserve
RUN R -e 'install.packages(c("Rserve"),repos="https://rforge.net/")'
# install R packages sources
RUN R -e 'install.packages(c("BiocManager", "remotes"),repos="http://cran.rstudio.com/")'
# Install R packages
RUN R -e "install.packages(c('shinythemes', 'shinymanager', 'here', 'ggplot2', 'dplyr', 'tidyr', 'stringr','gprofiler2', 'plotly', 'rentrez', 'HGNChelper', 'DT'), dependencies=TRUE)"
RUN R -e 'BiocManager::install(c("mygene"))'

# Set the working directory in the container
WORKDIR /app

# Copy the application files into the container
COPY app.R /app/
  COPY data-preparation.R /app/
  
  # Run the Shiny app
  CMD ["R", "-e", "shiny::runApp('/app/app.R', host='0.0.0.0', port=3838)"]

# Expose the port that the app will run on
EXPOSE 3838
