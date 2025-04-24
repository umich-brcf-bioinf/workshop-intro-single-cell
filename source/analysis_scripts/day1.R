# Create project directories
dir.create('scripts', showWarnings = FALSE, recursive = TRUE)
dir.create('results/figures', showWarnings = FALSE, recursive = TRUE)
dir.create('results/tables', showWarnings = FALSE, recursive = TRUE)
dir.create('results/rdata', showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)

# Modifies the size an R object can be in the session
options(future.globals.maxSize = 1e9)

# -------------------------------------------------------------------------
# To load data from 10X Cell Ranger

# Collect the input directories
#   Each sample dir contains barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#   Naming the sample_dirs vector makes Seurat name the
#     samples in the corresponding manner, which is nice for us.
sample_dirs = list(
  HODay0replicate1  = "inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate1", 
  HODay0replicate2  = "inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate2",
  HODay0replicate3  = "inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate3",
  HODay0replicate4  = "inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate4",
  HODay7replicate1  = "inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate1",
  HODay7replicate2  = "inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate2",
  HODay7replicate3  = "inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate3",
  HODay7replicate4  = "inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate4",
  HODay21replicate1 = "inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate1",
  HODay21replicate2 = "inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate2",
  HODay21replicate3 = "inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate3",
  HODay21replicate4 = "inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate4")

# Create the expression matrix from sample dirs
#   Read10X needs a *vector* instead of a *list*, so we use *unlist* to convert
geo_mat = Read10X(data.dir = unlist(sample_dirs))

# -------------------------------------------------------------------------
# Build BPCells input dir
# To use BPCells (to save some memory), you can transform 
# the expression matrix data structure above into BPCells files.
# BPCells uses these files for processing, but you typically never look at 
# their contents directly
write_matrix_dir(mat = geo_mat, dir = '~/ISC_R/bpcells', overwrite = TRUE)

# -------------------------------------------------------------------------
# Cleanup
# WORKSHOP SPECIFIC
# Since we'll now be reading in from BPCells files, we will remove geo_mat
# from the environment and then prompt RStudio to run a "garbage collection"
# to free up unused memory
rm(geo_mat)
gc()

# -------------------------------------------------------------------------
# Create expression matrix and Seurat object from BPCells files
geo_mat = open_matrix_dir(dir = '~/ISC_R/bpcells')

# -------------------------------------------------------------------------
# Create seurat object
geo_so = CreateSeuratObject(counts = geo_mat, min.cells = 1, min.features = 50)
geo_so

# -------------------------------------------------------------------------
# Examine Seurat object
head(geo_so@meta.data)

# -------------------------------------------------------------------------
# Save the Seurat object
saveRDS(geo_so, file = 'results/rdata/geo_so_unfiltered.rds')

# Practice shutting down session and reloading the R object
# WORKSHOP SPECIFIC
# -------------------------------------------------------------------------
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
options(future.globals.maxSize = 1e9)

# -------------------------------------------------------------------------
# Load the Seurat object
geo_so = readRDS('results/rdata/geo_so_unfiltered.rds')


# -------------------------------------------------------------------------
# Read sample attributes
# Load the expanded phenotype columns (condition, day, and replicate)
phenos = read.csv('inputs/prepared_data/phenos.csv')
head(phenos, 3)

# =========================================================================
# Secondary QC and filtering
# =========================================================================

# -------------------------------------------------------------------------
# Examine Seurat metdata

# We could do this to show metadata in the console pane:
head(geo_so@meta.data)

# For readability, we'll use View() to open the result in a scrollable tab
View(head(geo_so@meta.data))

# -------------------------------------------------------------------------
# Make a temp table that joins Seurat loaded metadata with expanded phenotype columns
# Preserve the rownames and order samples based on order in phenos.csv.
tmp_meta = geo_so@meta.data %>% 
  rownames_to_column('tmp_rowname') %>%
  left_join(phenos, by = 'orig.ident') %>%
  mutate(orig.ident = factor(orig.ident, phenos$orig.ident)) %>%
  mutate(time = factor(time, unique(phenos$time))) %>%
  column_to_rownames('tmp_rowname')

View(head(tmp_meta))

# -------------------------------------------------------------------------
# Assign tmp_meta back to geo_so@meta.data and 
# reset the default identity cell name
geo_so@meta.data = tmp_meta
Idents(geo_so) = 'orig.ident'
head(geo_so@meta.data)

# -------------------------------------------------------------------------
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
ggsave(filename = 'results/figures/qc_nFeature_violin.png', width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Review count violin plots
VlnPlot(geo_so, features = 'nCount_RNA', assay = 'RNA', layer = 'counts') + 
  NoLegend()
ggsave(filename = 'results/figures/qc_nCount_violin.png', width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Consider mitochondrial transcripts
# We use "mt" because this is mouse, depending on the organism, this might need to be changed
geo_so$percent.mt = PercentageFeatureSet(geo_so, pattern = '^mt-')

# Use summary() to quickly check the range of values
summary(geo_so$percent.mt)

head(geo_so)

# -------------------------------------------------------------------------
# Review mitochondrial violin plots
VlnPlot(geo_so, features = 'percent.mt', assay = 'RNA', layer = 'counts') + 
  NoLegend() + 
  geom_hline(yintercept = 25) + 
  geom_hline(yintercept = 20) + 
  geom_hline(yintercept = 15) + 
  geom_hline(yintercept = 10)
ggsave(filename = 'results/figures/qc_mito_violin.png', width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Filter to exclude suspect cells and assign to new Seurat object
geo_so = subset(geo_so, subset = nFeature_RNA > 300 & percent.mt < 15)
geo_so

# -------------------------------------------------------------------------
# Examine remaining cell counts
cell_counts_post_tbl = geo_so@meta.data %>% count(orig.ident, name = 'postfilter_cells')
cell_counts_post_tbl

# -------------------------------------------------------------------------
# Show cell counts before and after filtering
cell_counts_tbl = cell_counts_pre_tbl %>% left_join(cell_counts_post_tbl, by = 'orig.ident')
cell_counts_tbl

# -------------------------------------------------------------------------
write_csv(cell_counts_tbl, file = 'results/tables/cell_filtering_counts.csv')

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