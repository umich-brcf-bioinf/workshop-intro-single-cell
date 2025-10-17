# =========================================================================
# Independent Exercise - Day 3 Startup
# =========================================================================

# After restarting our R session, load the required libraries
library(Seurat)
library(BPCells)
library(tidyverse)
library(scCATCH)

# Load in seurat object with alternative clustering results from yesterday's exercises
exso3 = readRDS('results/rdata/geo_so_sct_clustered_exercise.rds')
exso3 # check that object loaded

## NOTE - BEFORE STOPPING WORK ON THE EXERCISES REMEMBER TO POWER DOWN AND RESTART R SESSION !!!!


### Day 3 Exercise 1 - examine marker genes and cell type predictions for alternative clustering
## ----------------------------------------------------------
# Check what identities are set
Idents(exso3) %>% head()

## use same values as previous exercise
pcs = 6
res = 0.2 

## Set identities to clustering for selected resolution
exso3 = SetIdent(exso3, value = paste0('int.sct.rpca.clusters', res))
Idents(exso3) %>% head()


## ----------------------------------------------------------
## Generate cluster markers to see how that changes with new parameters 
exso3 = PrepSCTFindMarkers(exso3, assay = "SCT") # NOTE - this step will take some time to run
exso3_markers = FindAllMarkers(exso3, only.pos = TRUE)
head(exso3)

# Create table of top 5 markers per cluster (using default ranking)
top_5 = exso3_markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% slice_head(n = 5)
head(top_5, n = 10) # look at results
write_csv(top_5, file = paste0('results/tables/top5_marker_genes_exercise.csv'))


#### Day 3 Exercise 2 - Generate predictions for alternative clustering 
## ----------------------------------------------------
## Next - run scCATCH predictions for alternative clustering results
exso3_catch = createscCATCH(data = exso3@assays$SCT@counts, cluster = as.character(Idents(exso3)))
catch_markers = exso3_markers %>% rename('logfc' = 'avg_log2FC')
exso3_catch@markergene = exso3_markers
exso3_catch@marker = cellmatch[cellmatch$species == 'Mouse' & cellmatch$tissue %in% c('Blood', 'Peripheral Blood', 'Muscle', 'Skeletal muscle', 'Epidermis', 'Skin'), ]
exso3_catch = findcelltype(exso3_catch)

# Check predictions
exso3_catch@celltype %>% select(cluster, cell_type, celltype_score)


#### Day 3 Exercise 3 - Add predictions to Seurat object
## ------------------------------------------------------
## Use predictions to label clusters and UMAP plot
catch_celltypes = exso3_catch@celltype %>% select(cluster, cell_type)
colnames(catch_celltypes)[2] = paste0('cell_type.',pcs,'PC.',res,'res')
new_metadata = exso3@meta.data %>% 
  left_join(catch_celltypes, 
            by = c('seurat_clusters' = 'cluster')) # using `seurat_clusters`, which will store the most recently generated cluster labels for each cell
rownames(new_metadata) = rownames(exso3@meta.data) #  We are implicitly relying on the same row order!
exso3@meta.data = new_metadata # Replace the meta.data
head(exso3@meta.data)


#### Day 3 Exercise 4 - Plot UMAP with new cluster labels
## ------------------------------------------------------
catch_umap_plot = DimPlot(exso3, group.by = paste0('cell_type.',pcs,'PC.',res,'res'), 
                          label = TRUE, reduction = paste0('umap.integrated.sct.rpca_alt', res))
catch_umap_plot


# Save the plot to file
# output to file, including the number of PCs and resolution used to generate results
ggsave(filename = paste0('results/figures/umap_int_catch-labeled_', 
                         pcs,'PC.',res,'res','.png'),
       plot = catch_umap_plot, 
       width = 8, height = 6, units = 'in')


## ----------------------------------------------------------
## Clean up session 
rm(list=names(which(unlist(eapply(.GlobalEnv, is.ggplot))))); 
rm(catch_celltypes, catch_markers, exso3_catch, exso3_markers, new_metadata, top_5); 
gc()

## (Optional) - Save copy of exso3
saveRDS(exso3, file = paste0('results/rdata/geo_so_sct_integrated_with_markers_exercise.rds'))
rm(exso3)
gc()

## BEFORE PROCEEDING TO THE NEXT SECTION or closing window - POWER DOWN AND RESTART R SESSION

