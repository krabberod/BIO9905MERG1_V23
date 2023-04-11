### Under construction

We will use R (version 4.0.5 or later) and Rstudio (version 1.4.1 or later) in this course. 

Everybody should download and install R (https://www.r-project.org/), Rstudio (https://www.rstudio.com/) and the required packages before the course starts.

It is preferable that you have all packages installed before the course starts. For those of you having trouble with installing the packages, we will provide some help Monday afternoon, 4pm, but we will not be able to provide comprehensive IT-service though.

These are the required packages for the hands-on exercises of the course: 
```
install.packages("readxl")    # To read excel files
install.packages("tidyverse") # To manipulate and visualize data
install.packages("metacoder")
install.packages("magrittr") 

if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq","Biostrings","PCAtools"))
```
More packages, including those that have been mentione, but not necessarily demonstrated can be found in the script [R packages](Install_packages.R).