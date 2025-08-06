3 + 2

2^3

# =========================================================================
# Getting Started with Seurat
# =========================================================================

# -------------------------------------------------------------------------
# Get current working directory
getwd()

# Set working directory to ISC_R
setwd('~/ISC_R')

# -------------------------------------------------------------------------
# Create directory structure

dir.create('scripts', recursive = TRUE, showWarnings = FALSE)
dir.create('results/figures', recursive = TRUE, showWarnings = FALSE)
dir.create('results/tables', recursive = TRUE, showWarnings = FALSE)
dir.create('results/rdata', recursive = TRUE, showWarnings = FALSE)

# Load the necessary libraries
library(Seurat)
library(BPCells)
library(tidyverse)

options(future.globals.maxSize = 1e9)

# -------------------------------------------------------------------------
# To load data from 10X Cell Ranger

# Collect the input directories
#   Each sample dir contains barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#   Naming the sample_dirs vector makes Seurat name the
#     samples in the corresponding manner, which is nice for us.
sample_dirs = list(
  HODay0replicate1  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate1", 
  HODay0replicate2  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate2",
  HODay0replicate3  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate3",
  HODay0replicate4  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate4",
  HODay7replicate1  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate1",
  HODay7replicate2  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate2",
  HODay7replicate3  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate3",
  HODay7replicate4  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate4",
  HODay21replicate1 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate1",
  HODay21replicate2 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate2",
  HODay21replicate3 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate3",
  HODay21replicate4 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate4")

# Create the expression matrix from sample dirs
#   Read10X needs a *vector* instead of a *list*, so we use *unlist* to convert
geo_mat = Read10X(data.dir = unlist(sample_dirs))

# Build the BPCells input directory and write the sparse matrix to disk
write_matrix_dir(mat = geo_mat, dir = '~/ISC_R/bpcells', overwrite = TRUE)

# Cleanup and garbage collection
rm(geo_mat)
gc()

# Read BPCells version of the count matrix
geo_mat = open_matrix_dir(dir = '~/ISC_R/bpcells')

# Create the Seurat object from geo_mat
geo_so = CreateSeuratObject(count = geo_mat, min.cells = 1, min.features = 50)
geo_so

# Examine the meta.data slot of the Seurat
head(geo_so@meta.data)

# Save the Seurat object
saveRDS(geo_so, file = 'results/rdata/geo_so_unfiltered.rds')

# We practiced powering down the session to clear our memory usage, now restart
# Let's load everything again
# -------------------------------------------------------------------------
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
options(future.globals.maxSize = 1e9)

setwd("~/ISC_R/")

# -------------------------------------------------------------------------
# Load the Seurat object
geo_so = readRDS('results/rdata/geo_so_unfiltered.rds')

# Peek at the first few rows of meta.data
head(geo_so@meta.data)

# Read in the prepared phenotype file and add it to meta.data
phenos = read.csv('inputs/prepared_data/phenos.csv')
phenos

# -------------------------------------------------------------------------
# Make a temp table that joins Seurat loaded metadata with expanded phenotype columns
# Preserve the rownames and order samples based on order in phenos.csv.
tmp_meta = geo_so@meta.data %>% 
  rownames_to_column('tmp_rowname') %>%
  left_join(phenos, by = 'orig.ident') %>%
  mutate(orig.ident = factor(orig.ident, phenos$orig.ident)) %>%
  mutate(time = factor(time, unique(phenos$time))) %>%
  column_to_rownames('tmp_rowname')

head(tmp_meta)

# Check that the number of rows is the same
dim(tmp_meta)
dim(geo_so@meta.data)

# Reassign meta.data to be tmp_meta
geo_so@meta.data = tmp_meta
Idents(geo_so) = 'orig.ident'

# Verify that the substitution worked
head(geo_so@meta.data)

# Let's take a look at the number of cells per sample
cell_counts_pre_tbl = geo_so@meta.data %>% count(orig.ident, name = 'prefilter_cells')
cell_counts_pre_tbl

# -------------------------------------------------------------------------
# Review feature violin plots
VlnPlot(geo_so, features = 'nFeature_RNA', assay = 'RNA', layer = 'counts') + 
  NoLegend() + 
  geom_hline(yintercept = 500) + 
  geom_hline(yintercept = 400) + 
  geom_hline(yintercept = 300) + 
  geom_hline(yintercept = 200)
ggsave(filename = 'results/figures/qc_nFeature_violin.png',
       width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Review count violin plots
VlnPlot(geo_so, features = 'nCount_RNA', assay = 'RNA', layer = 'counts') + 
  NoLegend()
ggsave(filename = 'results/figures/qc_nCount_violin.png',
       width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Consider mitochondrial transcripts
# We use "mt" because this is mouse, depending on the organism, this might need to be changed
geo_so$percent.mt = PercentageFeatureSet(geo_so, pattern = '^mt-')

# Use summary() to quickly check the range of values
summary(geo_so$percent.mt)

# -------------------------------------------------------------------------
# Review mitochondrial violin plots
VlnPlot(geo_so, features = 'percent.mt', assay = 'RNA', layer = 'counts') + 
  NoLegend() + 
  geom_hline(yintercept = 25) + 
  geom_hline(yintercept = 20) + 
  geom_hline(yintercept = 15) + 
  geom_hline(yintercept = 10)
ggsave(filename = 'results/figures/qc_mito_violin.png',
       width = 12, height = 6, units = 'in')

# Based on the filtering in the paper how many cells would remain?
# nFeature_RNA > 500 and percent.mt < 25
geo_so@meta.data %>% filter(nFeature_RNA > 500 & percent.mt < 25) %>% count(orig.ident)

cell_counts_pre_tbl

# Let's subset our geo_so so that nFeature_RNA > 300 & percent.mt < 15
geo_so = subset(geo_so, subset = nFeature_RNA > 300 & percent.mt < 15)
geo_so

cell_counts_post_tbl = geo_so@meta.data %>% count(orig.ident, name = 'postfilter_cells')
cell_counts_post_tbl

# -------------------------------------------------------------------------
# Show cell counts before and after filtering
cell_counts_tbl = cell_counts_pre_tbl %>% left_join(cell_counts_post_tbl, by = 'orig.ident')
cell_counts_tbl

# -------------------------------------------------------------------------
# Save the current Seurat object
saveRDS(geo_so, file = 'results/rdata/geo_so_filtered.rds')

# =========================================================================
# Normalization
# =========================================================================

# -------------------------------------------------------------------------
# Separate sample data into layers
geo_so[['RNA']] = split(geo_so[['RNA']], f = geo_so$orig.ident)
geo_so

# -------------------------------------------------------------------------
# Normalize the data with SCTransform
geo_so = SCTransform(geo_so)
geo_so

# -------------------------------------------------------------------------
# Save Seurat object
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_normalized.rds')
