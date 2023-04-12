# BIO9905MERG1_V23
# Please make sure that you have the required R packages installed.
# You can get a list of the package already installed on your computer s by executing

installed.packages()

# in R alternatively you can just run the installation for each package
# to make sure that you have the latest version.
# Run these installation commands line-by-line in R (or Rstudio)
# and answer yes if you are asked to update any previously installed packages:


install.packages("readxl")    # To read excel files
install.packages("metacoder") 
install.packages("tidyverse") # To manipulate and visualize data
# Tidyverse is a collection of packages that are really useful for data manipulation,
# an visualization. These packages are included and you do not need to install them 
# separately (although you can if you want to)
# install.packages("dplyr")     # To manipulate dataframes
# install.packages("tibble")    # To work with data frames
# install.packages("tidyr")     # To work with data frames
# install.packages("stringr")   # To manipulate strings
# install.packages("ggplot2")   # To do plots
# install.packages("readr")     # To read and write files
# install.packages("kableExtra")

install.packages("magrittr")

# Some packages, (for instance dada2 and phyloseq) are in a different repository
# called Bioconductor. To install package from the Bioconductor repository
# you need to use the BiocManager package: 
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq","Biostrings","PCAtools"))

install.packages("devtools")               # Developer tools
devtools::install_github("tobiasgf/lulu")

# If you run into problems when running the above line try to install one package
# at the time. i.e.: 
# if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
# BiocManager::install("dada2")

# install.packages("Matrix")	# Manipulation of large matrices
# install.packages("wTO") 	# Newtork analysis
# install.packages("igraph")	# Network analysis
# install.packages("standardize")
# install.packages("vegan")

# Packages related to ecological analysis. Only shown as examples.
# install.packages("spaa")           # Installs the ecological package spaa
# install.packages("compositions")   # To work with compositional data
# install.packages("zCompositions")  # To work with compositional data
# install.packages("mixOmics")       # Multivariate methods
# install.packages("ape")            # Phylogenetic tools
# install.packages("recluster")      # Clustering tools
# install.packages("dendextend")     # To work with dendrograms
# install.packages("corrplot")       # makes nice correlation plots
# install.packages("RcmdrMisc")      # diverse tools

# BiocManager::install("microbiome")


# Network packages
# library(devtools)
# install_github("zdk123/SpiecEasi")
# BiocManager::install("SpiecEasi") # Network construction
# devtools::install_github("pr2database/pr2database") # Installs directly from github resources that are not in R repos
# devtools::install_github("GuillemSalazar/EcolUtils") # Installs other tools for ecological analyses
# devtools::install_github('fawda123/ggord')
# This package might cause a problem. Here are some possible solutions:
# https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/
# https://stackoverflow.com/questions/37776377/error-when-installing-an-r-package-from-github-could-not-find-build-tools-neces

