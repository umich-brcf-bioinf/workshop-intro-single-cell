# =========================================================================
# Day 2 set up
# =========================================================================

# -------------------------------------------------------------------------
# Set working directory
# Confirm you are in the right working directory; should be ~/ISC_R
getwd()
# Reset if necessary
setwd("~/ISC_R")

# -------------------------------------------------------------------------
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
options(future.globals.maxSize = 1e9)


# Load data  ---------------------------------------------------
geo_so = readRDS('~/ISC_R/inputs/prepared_data/rdata/geo_so_sct_normalized.rds')
geo_so

# Alternatively, this starts from the data you saved yesterday
#geo_so = readRDS('~/ISC_R/results/rdata/geo_so_sct_normalized.rds')

# =========================================================================
# PCA and Integration
# =========================================================================

# Build a PCA and add it to the Seurat object
geo_so = RunPCA(geo_so, reduction.name = 'unintegrated.sct.pca')
geo_so

# -------------------------------------------------------------------------
# look at first PC alone first
heatmap_1 <- DimHeatmap(geo_so, dims=1, cells=500, balanced=TRUE, fast = FALSE, reduction = 'unintegrated.sct.pca') # doesn't include which PC(s) are included in plot
heatmap_1 + 
  labs(title="PC1: top 30 genes x top 500 cells") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

# save to file
ggsave(filename = 'results/figures/qc_pca_heatmap-pc1.png',
       width = 6, height = 4.5, units = 'in')


#  -------------------------------------------------------
# plot cell by gene heatmaps for larger set of dimensions
DimHeatmap(geo_so, dims=1:18, cells=500, balanced=TRUE, combine = FALSE, reduction = 'unintegrated.sct.pca')

# need to use png() to write to file because this isn't a ggplot
png(filename = 'results/figures/qc_pca_heatmap-pc1-18.png', width = 12, height = 30, units = 'in', res = 300)
DimHeatmap(geo_so, dims=1:18, cells=500, balanced=TRUE, reduction = 'unintegrated.sct.pca')
dev.off()

# ------------------------------------------------------------
# Plot cell by gene heatmaps for a subset of dimensions 
DimHeatmap(geo_so, dims=c(1,5,10,20,30,50), cells=500, balanced=TRUE, reduction = 'unintegrated.sct.pca', nfeatures=20)

png(filename = 'results/figures/qc_pca_heatmap-select.png', width = 12, height = 10, units = 'in', res = 300)
DimHeatmap(geo_so, dims=c(1,5,10,20,30,50), cells=500, balanced=TRUE, reduction = 'unintegrated.sct.pca', nfeatures=20)
dev.off()

# Use the ? to see the documentation page for a function
?DimPlot

# -------------------------------------------------------------------------
# Visualize PCA (PC1 and PC2, unintegrated)
# first label by sample
DimPlot(geo_so, reduction = 'unintegrated.sct.pca', group.by = 'orig.ident') 

# then label by day, which looks more useful
head(geo_so@meta.data)
DimPlot(geo_so, reduction = 'unintegrated.sct.pca', group.by = 'time')

# save plot labeled by day
ggsave(filename = 'results/figures/qc_pca_plot_unintegrated_sct_day.png',
       width = 7, height = 6, units = 'in')

# Check documentation for integrate layers command
?IntegrateLayers

# -------------------------------------------------------------------------
# Load integrated data from prepared file
geo_so = readRDS('inputs/prepared_data/rdata/geo_so_sct_integrated.rds')
geo_so

# -------------------------------------------------------------------------
# Visualize PCA (PC1 and PC2, integrated)
DimPlot(geo_so, reduction = 'integrated.sct.rpca', group.by = 'time')
ggsave(filename = 'results/figures/qc_pca_plot_integrated_sct_day.png', 
       width = 7, height = 6, units = 'in')


# -------------------------------------------------------------------------
# Remove plot variables from the environment to avoid excessive memory usage
plots = c("heatmap_1")

# Only remove plots that actually exist in the environment
rm(list=Filter(exists, plots))
gc()

# -------------------------------------------------------------------------
# Save updated seurat object to file
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated.rds')

## Clustering and Projection

# -------------------------------------------------------------------------
# Visualize how many PCs to include using an elbow plot
ElbowPlot(geo_so, ndims = 50, reduction = 'unintegrated.sct.pca')
ggsave(filename = 'results/figures/qc_sct_elbow_plot.png',
       width = 8, height = 8, units = 'in')

# -------------------------------------------------------------------------
# Estimate the number of PCs to use for clustering with a function
check_pcs = function(so, reduction) {
  # quantitative check for number of PCs to include
  pct = so@reductions[[reduction]]@stdev / sum(so@reductions[[reduction]]@stdev) * 100
  cum = cumsum(pct)
  co1 = which(cum > 90 & pct < 5)[1]
  co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs = min(co1, co2) 
  
  return(pcs)
}

# Apply function to our data
pcs = check_pcs(geo_so, 'unintegrated.sct.pca')
pcs

# -------------------------------------------------------------------------
# Based some behind the scenes testing, we'll modify the number of PCs to use for clustering
pcs = 10

# =========================================================================
# Clustering and Projection
# =========================================================================

# Create KNN graph for PCs with `FindNeighbors()`
geo_so = FindNeighbors(geo_so, dims = 1:pcs, reduction = 'integrated.sct.rpca')

# Then generate clusters
geo_so = FindClusters(geo_so,
                      resolution = 0.4,
                      cluster.name = 'integrated.sct.rpca.clusters')

# look at meta.data to see cluster labels
head(geo_so@meta.data)


# -------------------------------------------------------------------------
# Create UMAP reduction
geo_so = RunUMAP(geo_so,
                 dims = 1:pcs,
                 reduction = 'integrated.sct.rpca',
                 reduction.name = 'umap.integrated.sct.rpca')

# Note a third reduction has been added: `umap.integrated.sct.rpca`
geo_so 

# Check parameters for `dim`
1:pcs

# -------------------------------------------------------------------------
# Visualize UMAP clusters: ID labels
post_integration_umap_plot_clusters = 
  DimPlot(geo_so, 
          group.by = 'seurat_clusters',
          label = TRUE,
          reduction = 'umap.integrated.sct.rpca') + 
  NoLegend()
post_integration_umap_plot_clusters

ggsave(filename = 'results/figures/umap_integrated_sct_clusters.png',
       plot = post_integration_umap_plot_clusters,
       width = 6, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Visualize UMAP clusters: clusters with labels, split by condition
post_integration_umap_plot_split_clusters = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'time',
          label = TRUE,
          reduction = 'umap.integrated.sct.rpca') + 
  NoLegend()
post_integration_umap_plot_split_clusters

ggsave(filename = 'results/figures/umap_integrated_sct_split_clusters.png',
       plot = post_integration_umap_plot_split_clusters, 
       width = 14, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Visualize UMAP clusters: day labels
post_integration_umap_plot_day = 
  DimPlot(geo_so,
          group.by = 'time',
          label = FALSE,
          combine = FALSE,
          reduction = 'umap.integrated.sct.rpca')[[1]] + 
  ggtitle("Post-integration clustering, by timepoint")
post_integration_umap_plot_day

ggsave(filename = 'results/figures/umap_integrated_sct_day.png',
       plot = post_integration_umap_plot_day,
       width = 8, height = 6, units = 'in')


# -------------------------------------------------------------------------
# Make tables of cell counts
clusters_counts_condition <- geo_so@meta.data %>% 
  select(time, integrated.sct.rpca.clusters) %>%
  summarise(counts = n(), .by = c(time, integrated.sct.rpca.clusters)) %>% 
  pivot_wider(names_from = integrated.sct.rpca.clusters, values_from = counts, names_sort=TRUE) %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "TOTAL"))
clusters_counts_condition 

# cells per cluster per sample
clusters_counts_sample <- geo_so@meta.data %>% 
  select(orig.ident, integrated.sct.rpca.clusters) %>%
  summarise(counts = n(), .by = c(orig.ident, integrated.sct.rpca.clusters)) %>% 
  pivot_wider(names_from = integrated.sct.rpca.clusters, values_from = counts, names_sort=TRUE) %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(., na.rm = TRUE) else "TOTAL"))
clusters_counts_sample

# -------------------------------------------------------------------------
# Remove plot variables from the environment to manage memory usage
plots = c("pre_integration_umap_plot_day", 
          "post_integration_umap_plot_clusters", 
          "post_integration_umap_plot_split_clusters", 
          "post_integration_umap_plot_day")

# Only remove plots that actually exist in the environment
rm(list=Filter(exists, plots))
gc()

# -------------------------------------------------------------------------
# Save Seurat object
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_clustered.rds')
