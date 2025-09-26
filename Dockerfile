# Base image with R
FROM rocker/r-ver:4.3.2

# System dependencies for rstan and general R packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    gfortran \
    libblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up rstan options for compilation
ENV RSTAN_OPTIONS="auto_write = TRUE"
ENV CXX14FLAGS="-O3 -march=native -mtune=native"

# Install remotes package to install your GitHub package
RUN R -e "install.packages('remotes')"

# Install package dependencies
RUN R -e "install.packages(c(    
    'data.table',
    'DirichletReg',
    'dplyr',
    'ggh4x',
    'ggplot2',
    'grid',
    'methods',
    'rlang',
    'rstan',
    'stats',
    'stringr',
    'tidyr',
    'utils'))"

# Copy your package into the container
COPY . /townlet

# Set working directory
WORKDIR /townlet

# Install your package
RUN R -e "remotes::install_local('.')"

# Default command
CMD ["R"]
