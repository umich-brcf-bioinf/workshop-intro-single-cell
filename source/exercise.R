# Exercise: adjust number of PCs and resolution to see impact on clustering

# First - clear current Seurat object to free up memory & remove current results
rm(geo_so) 

# Then  load in integrated data & reset object on each iteration to avoid exceeding allocated space
geo2_so = readRDS('/home/workshop/rcavalca/ISC_R/results/rdata/geo_so_sct_integrated.rds')

# look at elbow plot
ElbowPlot(geo2_so, ndims = 50, reduction = 'unintegrated.sct.pca')

## Clustering
# Round 1: manually adjust number of PCs to include in clustering
#pcs = 20 # increase number of PCs
pcs = 10 # reduce number of PCs

# generate nearest neighbor (graph), using selected number of PCs
geo2_so = FindNeighbors(geo2_so, dims = 1:pcs, reduction = 'integrated.sct.rpca')

# Round 2: adjust resolution after testing PCs (remember this only impacts the how the boundaries of the neighbors are drawn, not the underlying NN graph/structure)
res = 0.4
# res = 1.0
# res = 0.2

# generate clusters, using 
geo2_so = FindClusters(geo2_so, resolution = res, cluster.name = 'integrated.sct.rpca.clusters')

# look at meta.data to see cluster labels
head(geo2_so@meta.data)

# Prep for UMAP
geo2_so = RunUMAP(geo2_so, dims = 1:pcs, reduction = 'integrated.sct.rpca', 
                 reduction.name = 'umap.integrated.sct.rpca')
geo2_so

# look at clustering results
post_integration_umap_clusters_testing = 
  DimPlot(geo2_so, group.by = 'seurat_clusters', label = TRUE, 
          reduction = 'umap.integrated.sct.rpca') + NoLegend()
post_integration_umap_clusters_testing # look at plot

# output to file, including the number of PCs and resolution used to generate results
ggsave(filename = paste0('results/figures/umap_integrated_sct_clusters', 
                         pcs,'PCs',res,'res.png'),
       plot = post_integration_umap_plot_clusters, 
       width = 8, height = 6, units = 'in')

## generate markers and annotate clusters to see how that changes

# prep for cluster comparisons
geo2_so = SetIdent(geo2_so, value = 'integrated.sct.rpca.clusters')
geo2_so = PrepSCTFindMarkers(geo2_so)

# run comparisons for each cluster to generate markers
geo2_markers = FindAllMarkers(geo2_so, only.pos = TRUE)


# run cell type predictions for current clustering results
library(scCATCH)

# create scCATCH object, using count data
geo2_catch = createscCATCH(data = geo2_so@assays$SCT@counts, cluster = as.character(Idents(geo2_so)))

# add marker genes to use for predictions
catch_markers = geo2_markers %>% rename('logfc' = 'avg_log2FC')
geo2_catch@markergene = geo2_markers

# specify tissues/cell-types from the scCATCH reference
geo2_catch@marker = cellmatch[cellmatch$species == 'Mouse' & cellmatch$tissue %in% c('Blood', 'Peripheral Blood', 'Muscle', 'Skeletal muscle', 'Epidermis', 'Skin'), ]

# run scCATCH to generate predictions
geo2_catch = findcelltype(geo2_catch)

# look at the predictions
geo2_catch@celltype %>% select(cluster, cell_type, celltype_score)

## annotate clusters
# Extract the cell types only to merge into the meta.data
catch_celltypes = geo2_catch@celltype %>% select(cluster, cell_type)

# Merge cell types in but as a new table to slide into @meta.data
new_metadata = geo2_so@meta.data %>% left_join(catch_celltypes, 
                                              by = c('integrated.sct.rpca.clusters' = 'cluster'))
rownames(new_metadata) = rownames(geo2_so@meta.data) #  We are implicitly relying on the same row order!

# Replace the meta.data
geo2_so@meta.data = new_metadata 
head(geo2_so@meta.data)

catch_umap_plot = DimPlot(geo2_so, group.by = 'cell_type', 
                          label = TRUE, reduction = 'umap.integrated.sct.rpca')
catch_umap_plot



#########

## Extension - how might you generate DE comparisons between D21 and D7 for all annotated clusters?

