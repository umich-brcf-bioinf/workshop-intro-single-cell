# =========================================================================
# Independent Exercise - Day 2 Startup
# =========================================================================

# After restarting our R session, load the required libraries & input data
library(Seurat)
library(BPCells)
library(tidyverse)

# Use provided copy of integrated data
exso2 = readRDS('inputs/prepared_data/rdata/geo_so_sct_integrated.rds')
exso2 # check that object loaded


### Day 2 Exercise 1 - Clustering with reduced number of PCs
# -------------------------------------------------------------------------
# Testing fewer PCs for clustering 

# look at elbow plot to check PCs and consider alternatives to the number selected for main materials
ElbowPlot(exso2, ndims = 50, reduction = 'unintegrated.sct.pca')

# select alternative value to try (can choose a number <10 or >10 PCs)
pcs = 6

# generate nearest neighbor (graph), using selected number of PCs
exso2 = FindNeighbors(exso2, dims = 1:pcs, reduction = 'integrated.sct.rpca')

# -------------------------------------------------------------------------
# start with one resolution
res = 0.2

# generate clusters, using `pcs` and `res` to make a custom cluster name that will be added to the metadata
exso2 = FindClusters(exso2, resolution = res, 
                     cluster.name = paste0('int.sct.rpca.clusters', res))

# look at meta.data to see cluster labels
head(exso2@meta.data)

# run UMAP reduction to prepare for plotting
exso2 = RunUMAP(exso2, dims = 1:pcs, 
                reduction = 'integrated.sct.rpca', 
                reduction.name = paste0('umap.integrated.sct.rpca_alt', res))

# check object to see if named reduction was added
exso2


## Challenge 1 - solution option
# -------------------------------------------------------------------------
# use for loop to generate clustering with alternative resolutions
for(i in c(0.4, 0.8, 1.0)){
  exso2 = FindClusters(exso2, resolution = i, 
                       cluster.name = paste0('int.sct.rpca.clusters', i))
  
  # look at meta.data to see cluster labels
  head(exso2@meta.data)
  
  # run UMAP reduction to prepare for plotting
  exso2 = RunUMAP(exso2, dims = 1:pcs, 
                  reduction = 'integrated.sct.rpca', 
                  reduction.name = paste0('umap.integrated.sct.rpca_alt', i))
  
  # check object to see if multiple cluster resolutions are added
  head(exso2@meta.data)
}



### Day 2 Exercise 2 - Plotting alternative clustering results
# -------------------------------------------------------------------------
# plot clustering results
post_integration_umap_clusters_testing = 
  DimPlot(exso2, group.by = paste0('int.sct.rpca.clusters', res), label = TRUE, 
          reduction = paste0('umap.integrated.sct.rpca_alt', res)) + NoLegend()
post_integration_umap_clusters_testing # look at plot

# output to file, including the number of PCs and resolution used to generate results
ggsave(filename = paste0('results/figures/umap_int_sct_clusters_exercise_', 
                         pcs,'PC.',res,'res','.png'),
       plot = post_integration_umap_clusters_testing, 
       width = 8, height = 6, units = 'in')


## Challenge 2 - solution option
# -------------------------------------------------------------------------
# use for loop to visualize clustering across tested resolutions
post_integration_umap_plots <- c()
for(i in c(0.4, 0.8, 1.0)){
  res_type = paste0("res_", i)
  post_integration_umap_plots[[res_type]] = 
    DimPlot(exso2, group.by = paste0('int.sct.rpca.clusters', i), label = TRUE, 
            reduction = paste0('umap.integrated.sct.rpca_alt', i)) + NoLegend()
}

# look at plots for each resolution stored in list
post_integration_umap_plots
# remove plot list to clean up session 
rm(post_integration_umap_plots)


## ----------------------------------------------------------
## Clean up session, including any plot objects
rm(list=names(which(unlist(eapply(.GlobalEnv, is.ggplot)))));
gc()

## Save copy of Seurat object in current state to file
saveRDS(exso2, file = paste0('results/rdata/geo_so_sct_clustered_exercise.rds'))
rm(exso2)
gc()

