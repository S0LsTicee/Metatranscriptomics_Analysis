#This code detects both sample and taxon level outliers using PCA-based (Mahalanobis distance) method. 
# Install Packages
my_packages <- c("tidyverse", "vegan", "ggplot2", "ggrepel",
                   "RColorBrewer", "patchwork", "factoextra", "mvoutlier")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(factoextra)
library(mvoutlier)

#Load Files
