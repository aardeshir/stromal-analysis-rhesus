# Stromal Cell Analysis in Primate Colon
# Author: A. Ardeshir (Ardeshirlab - aardeshir@tulane.edu)
# Date: 2025-01-08

# Setup and Data Loading
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(scater)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(gridExtra)
})

# Create output directories
dir.create(file.path('output', 'stromal_analysis', 'plots'), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path('output', 'stromal_analysis', 'processed'), 
           recursive = TRUE, showWarnings = FALSE)

# Load data
metadata_add <- read.csv(file.path('data', 'Additional_metadata.csv'), header = TRUE)
gcm_PR03 <- readRDS(file.path('PR03_Chronic', 'data', 'processed', 'counts.RDS'))
locs_PR03 <- readRDS(file.path('PR03_Chronic', 'data', 'processed', 'xy.RDS'))
negcounts_PRO3 <- readRDS(file.path('PR03_Chronic', 'data', 'processed', 'negcounts.RDS'))
metadata_PR03 <- readRDS(file.path('PR03_Chronic', 'data', 'processed', 'metadata.RDS'))