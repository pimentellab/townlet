# Base image: R 4.3.2 with RStudio Server
FROM rocker/rstudio:4.3.2

# Install system dependencies for rstan and general R packages
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure RStan compilation flags
RUN mkdir -p /home/rstudio/.R && \
    echo "CXX14=g++" >> /home/rstudio/.R/Makevars && \
    echo "CXX14FLAGS=-O3 -march=native -mtune=native" >> /home/rstudio/.R/Makevars && \
    echo "CXX11FLAGS=-O3 -march=native -mtune=native" >> /home/rstudio/.R/Makevars && \
    chown -R rstudio:rstudio /home/rstudio/.R

# Install remotes and your package (with dependencies)
RUN R -e "install.packages('remotes')" \
    && R -e "remotes::install_github('pimentellab/townlet')"

# Expose RStudio Server default port
EXPOSE 8787

# Default user is already set to 'rstudio' with password 'rstudio'
# (you can override with environment variables at run time)
