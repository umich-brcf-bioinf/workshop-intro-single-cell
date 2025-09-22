# =========================================================================
# Independent Exercise - Day 1 Startup
# =========================================================================

# After restarting our R session, load the required libraries & input data
library(Seurat)
library(BPCells)
library(tidyverse)

setwd('~/ISC_R')

# Load the unfiltered version and give it a new variable name
exso = readRDS('./inputs/prepared_data/rdata/geo_so_unfiltered.rds')

# Add percent.mt column to meta.data
exso$percent.mt = PercentageFeatureSet(exso, pattern = '^mt-')

## NOTE - BEFORE STOPPING WORK ON THE EXERCISES REMEMBER TO POWER DOWN AND RESTART R SESSION !!!!

# Exercise 1
exso = subset(exso, nFeature_RNA >= 1000)

ex1 = exso@meta.data %>% count(orig.ident, name = 'postfilter_cells')
ex1

# Exercise 2
exso = subset(exso, nCount_RNA >= 1000)

ex2 = exso@meta.data %>% count(orig.ident, name = 'postfilter_cells')
ex2

# Exercise 3
exso = subset(exso, percent.mt < 10)

ex3 = exso@meta.data %>% count(orig.ident, name = 'postfilter_cells')
ex3

# Exercise 4
exso = subset(exso, nFeature_RNA >= 1000 & nCount_RNA >= 1000 & percent.mt < 10)

ex4 = exso@meta.data %>% count(orig.ident, name = 'postfilter_cells')
ex4
