##### Day 1 - Getting Started with Seurat

# Create project directories ----------------------------------------------
dir.create('scripts', showWarnings = FALSE, recursive = TRUE)
dir.create('results/figures', showWarnings = FALSE, recursive = TRUE)
dir.create('results/tables', showWarnings = FALSE, recursive = TRUE)
dir.create('results/rdata', showWarnings = FALSE, recursive = TRUE)

# Load libraries  ----------------------------------------------
library(Seurat)
library(BPCells)
library(tidyverse)

options(future.globals.maxSize = 1e9)

# Puts the data "on disk" rather than "in memory"  -------------
geo_mat = open_matrix_dir(dir = 'inputs/bpcells')

# Create seurat object  ----------------------------------------------
geo_so = CreateSeuratObject(counts = geo_mat, min.cells = 1, min.features = 50)
geo_so

# Examine Seurat object ----------------------------------------------
head(geo_so@meta.data)
  
# Save the Seurat object ----------------------------------------------
saveRDS(geo_so, file = 'results/rdata/geo_so_unfiltered.rds')

##### Day 1 - Secondary QC and filtering

# Examine Seurat metdata  ----------------------------------------------
head(geo_so@meta.data)

# Make metadata more granular ------------------------------------------
# Load the expanded phenotype columns (condition, day, and replicate)
phenos = read.csv('inputs/phenos.csv')
head(phenos, 3)

# Make a temp table that joins Seurat loaded metadata with expanded phenotype columns
# (Also, preserve the rownames and order samples based on order in phenos.csv)
tmp_meta = geo_so@meta.data %>% 
  rownames_to_column('tmp_rowname') %>%
  left_join(phenos, by = 'orig.ident') %>%
  mutate(orig.ident = factor(orig.ident, phenos$orig.ident)) %>%
  column_to_rownames('tmp_rowname')
head(tmp_meta)  

# Assign tmp_meta back to geo_so@meta.data and reset the default identity cell name
geo_so@meta.data = tmp_meta
Idents(geo_so) = 'orig.ident'

# Verify that the change was made with
head(geo_so@meta.data)

# Assign overall cell counts  ------------------------------------------
cell_counts_pre_tbl = geo_so@meta.data %>% count(orig.ident, name = 'prefilter_cells')
cell_counts_pre_tbl

# Review feature violin plots   ----------------------------------------
VlnPlot(geo_so, features = 'nFeature_RNA', assay = 'RNA', layer = 'counts') + NoLegend() + geom_hline(yintercept = 500) + geom_hline(yintercept = 400) + geom_hline(yintercept = 300) + geom_hline(yintercept = 200)
ggsave(filename = 'results/figures/qc_nFeature_violin.png', width = 12, height = 6, units = 'in')

# Review count violin plots   ------------------------------------------
VlnPlot(geo_so, features = 'nCount_RNA', assay = 'RNA', layer = 'counts') + NoLegend()
ggsave(filename = 'results/figures/qc_nCount_violin.png', width = 12, height = 6, units = 'in')

# Consider mitochondrial transcripts  ----------------------------------
# We use "mt" because this is mouse, depending on the organism, this might need to be changed
geo_so$percent.mt = PercentageFeatureSet(geo_so, pattern = '^mt-')

# Use summary() to quickly check the range of values
summary(geo_so$percent.mt)

# Preview the new meta.data table
head(geo_so@meta.data)

# Review mitochondrial violin plots   ----------------------------------
VlnPlot(geo_so, features = 'percent.mt', assay = 'RNA', layer = 'counts') + NoLegend() + geom_hline(yintercept = 25) + geom_hline(yintercept = 20) + geom_hline(yintercept = 15) + geom_hline(yintercept = 10)
ggsave(filename = 'results/figures/qc_mito_violin.png', width = 12, height = 6, units = 'in')

# Filter to exclude suspect cells and assign to new Seurat object ------
geo_so = subset(geo_so, subset = nFeature_RNA > 300 & percent.mt < 15)
geo_so

# Examine remaining cell counts ----------------------------------------
cell_counts_post_tbl = geo_so@meta.data %>% count(orig.ident, name = 'postfilter_cells')
cell_counts_post_tbl

# Show cell counts before and after filtering --------------------------
cell_counts_tbl = cell_counts_pre_tbl %>% left_join(cell_counts_post_tbl, by = 'orig.ident')
cell_counts_tbl
write_csv(cell_counts_tbl, file = 'results/tables/cell_filtering_counts.csv')

# Save the current Seurat object ---------------------------------------
saveRDS(geo_so, file = 'results/rdata/geo_so_filtered.rds')

# Load the Seurat object ----------------------------------------------
geo_so = readRDS('results/rdata/geo_so_filtered.rds')
