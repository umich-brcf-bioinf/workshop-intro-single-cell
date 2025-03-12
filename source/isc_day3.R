### Day 3 ---------

# load the libraries
library(Seurat)
library(BPCells)
library(tidyverse)

# Load data --------
geo_so = readRDS('/home/workshop/damki/ISC_R/results/rdata/geo_so_sct_clustered.rds')

##### Day 3 - Marker identification and visualization 

# Find empirical markers  -------------------------------------------------
# Prep for cluster comparisons
geo_so = SetIdent(geo_so, value = 'integrated.sct.rpca.clusters')
geo_so = PrepSCTFindMarkers(geo_so)

# Run comparisons for each cluster to generate markers
geo_markers = FindAllMarkers(geo_so, only.pos = TRUE, min.pct = 0.05)

# Write out full cluster marker results to file
write_csv(geo_markers, file = 'results/tables/marker_genes_0.4res.csv')

# Take a look at the first few rows of the result
head(geo_markers)


# Identify marker genes for each cluster ----------------------------------
# Create table of top 5 markers per cluster (using default ranking)
top_5 = geo_markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% slice_head(n = 5)

# Look at results
head(top_5, n = 10)


# Visualize top marker genes as dot plot ----------------------------------
top_5_sct_dot_plot = DotPlot(geo_so, features = unique(top_5$gene)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = 'Top 5 Cluster Genes by FDR and avg_log2FC') + coord_flip()
top_5_sct_dot_plot


# Save dot plot of top marker genes ---------------------------------------
ggsave(filename = 'results/figures/markers_top_5_sct_dot_plot.png', 
       plot = top_5_sct_dot_plot, 
       width = 8, height = 18, units = 'in') 

# Check mitochondrial gene expression -------------------------------------
percent_mito_plot = FeaturePlot(geo_so, features='percent.mt')

# save to file
ggsave(filename = 'results/figures/percent_umap_mito_plot.png', 
       plot = percent_mito_plot, 
       width = 6, height = 6, units = 'in')

percent_mito_plot

# Remove the plot variables from the environment to avoid excessive memory usage
plots = c("top_5_sct_dot_plot", 
          "top_5_rna_dot_plot", 
          "percent_mito_plot")

# Only remove plots that actually exist in the environment
rm(list=Filter(exists, plots))
gc()

# Save Seurat object and gene marker data ---------------------------------
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated_with_markers.rds')
saveRDS(geo_markers, file = 'results/rdata/geo_markers.rds')


##### Day 3 - Cell Type Annotation

# Load scCATCH  -----------------------------------------------------------
library(scCATCH)

# check that cell identities are set to expected resolution 
all(Idents(geo_so) == geo_so$integrated.sct.rpca.clusters)

# Annotate clusters using scCATCH  ----------------------------------------

# create scCATCH object, using count data
geo_catch = createscCATCH(data = geo_so@assays$SCT@counts, cluster = as.character(Idents(geo_so)))

# add marker genes to use for predictions
geo_catch@markergene = geo_markers

# specify tissues/cell-types from the scCATCH reference
geo_catch@marker = cellmatch[cellmatch$species == 'Mouse' & cellmatch$tissue %in% c('Blood', 'Peripheral Blood', 'Muscle', 'Skeletal muscle', 'Epidermis', 'Skin'), ]

# run scCATCH to generate predictions
geo_catch = findcelltype(geo_catch)

# look at the predictions
geo_catch@celltype %>% select(cluster, cell_type, celltype_score)

# Create lists of immune cells and associated gene markers --------------
immune_markers = list()
immune_markers[['Inflam. Macrophage']] = c('Cd14', 'Cxcl2') # Cd14 a- monocyte/macrophage cells
immune_markers[['Platelet']] = c('Pf4')
immune_markers[['Mast cells']] = c('Gata2', 'Kit')
immune_markers[['NK cells']] = c('Nkg7', 'Klrd1')
immune_markers[['B-cell']] = c( 'Ly6d', 'Cd19', 'Cd79b', 'Ms4a1')
immune_markers[['T-cell']] = c( 'Cd3d','Cd3e','Cd3g') # also Thy1

# Plot other immune to assist with cluster identification ---------------
immune_markers_plot = DotPlot(geo_so, features = immune_markers, assay = 'SCT')  +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 45, vjust = 0.5))

# save to file
ggsave(filename = 'results/figures/immune_markers_sct_dot_plot.png', 
       plot = immune_markers_plot, width = 10, height = 5, units = 'in')

immune_markers_plot

# Create lists of other cells and associated gene markers ---------------------------------------
other_markers = list()
other_markers[['Pericyte']] = c('Acan','Sox9')
other_markers[['SMC']] = c('Acta2', 'Myh11') # SMC = mesenchymal smooth-muscle cell/mesenchymal lineage
other_markers[['Keratinocytes']] = c('Thy1', 'Dlk1') # fibro progenitors aso=Thy1
other_markers[['Myofibroblasts']] = c('Tmem100', 'Cd34', 'Ly6c1') # hematopoetic stem/activated fibroblast=Cd34
other_markers[['Fibroblast']] = c('Dpt', 'Fn1', 'Col3a1')  # activated fib = Fn1
other_markers[['Endothelial']] = c('Pecam1', 'Cd38') # from wound healing; Pecam1 also exp in endothelial
other_markers[['HSC']] = c('Ltb', 'Cd74') # less well defined/conflicting definitions
other_markers[['Erythroid']] = c('Hba-a1')


# Plot known cell-type markers ---------------------------------------
other_markers_dot_plot = DotPlot(geo_so, features = other_markers, assay = 'SCT') +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 45, vjust = 0.5))

# save to file
ggsave(filename = 'results/figures/other_markers_sct_dot_plot.png', 
       plot = other_markers_dot_plot, width = 12, height = 5, units = 'in')

other_markers_dot_plot

# look at gene expression on UMAP plot
FeaturePlot(geo_so, features='Col12a1')


# Annotate clusters using modified predictions ----------------------------
# First - Extract the cell types only from the predictions
celltype_annos = geo_catch@celltype %>% select(cluster, cell_type) %>% 
  mutate(cluster = factor(cluster, levels = c(0:22))) %>% arrange(cluster)
celltype_annos

# Update annotations, remembering that cluster 0 = row 1 in table ---------
celltype_annos$cell_type[c(6,15,16)] <- "Inflammatory macrophage" # resolve cluster 5, 14, 15
celltype_annos$cell_type[c(12)] <- "Macrophage"
celltype_annos$cell_type[c(5,8,10)] <- "Platelet" # clusters 4,7,9

celltype_annos$cell_type[c(1,22)] <- "Pericyte"
celltype_annos$cell_type[c(2,9)] <- "Fibroblast" # revise clusters 1,8 based on markers
celltype_annos$cell_type[c(7)] <- "Myofibroblast" # revise cluster 6
celltype_annos$cell_type[c(4,14,17)] <- "Hematopoietic stem cell" # based on markers but could further revise
celltype_annos$cell_type[c(11,20)] <- "Mesenchymal stem/stromal cell" # based on Acta2 signal; cluster 10, 19
celltype_annos$cell_type[c(19)] <- "Erythroid" # cluster 18
celltype_annos$cell_type[c(23)] <- "Mast"

celltype_annos$cell_type[c(18, 21)] <- "Unknown" # since such small populations, reset cluster 17 & 20 as unknown for now----

# Merge cell types in but as a new table to slide into @meta.data ----------
copy_metadata = geo_so@meta.data
new_metadata = copy_metadata %>% left_join(celltype_annos, by = c('integrated.sct.rpca.clusters' = 'cluster'))
rownames(new_metadata) = rownames(geo_so@meta.data) #  We are implicitly relying on the same row order!

# Replace the meta.data
geo_so@meta.data = new_metadata 

head(geo_so@meta.data)

# Make a labeled UMAP plot of clusters ------------------------------------
catch_umap_plot = DimPlot(geo_so, group.by = 'cell_type', label = TRUE, reduction = 'umap.integrated.sct.rpca')
catch_umap_plot

ggsave(filename = 'results/figures/umap_integrated_catch.png', plot = catch_umap_plot, width = 10, height = 8, units = 'in')

catch_umap_condition_plot = DimPlot(geo_so, group.by = 'cell_type', split.by = 'day', label = TRUE, reduction = 'umap.integrated.sct.rpca')

ggsave(filename = 'results/figures/umap_integrated_catch_byCondition.png', 
       plot = catch_umap_plot, width = 10, height = 8, units = 'in')

catch_umap_condition_plot

# Discard all ggplot objects currently in environment ---------------------
# Ok since we saved the plots as we went along
rm(list=names(which(unlist(eapply(.GlobalEnv, is.ggplot))))); 
gc()

# Save Seurat object and annotations --------------------------------------
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated_with_catch.rds')
saveRDS(geo_catch, file = 'results/rdata/geo_catch.rds')

##### Day 3  - Differential Expression Analysis

# Make a day-celltype contrast --------------------------------------------

# set up combined label of day + celltype & assign as identities
geo_so$day.celltype = paste(geo_so$day, geo_so$cell_type, sep = '_')
# check labels
unique(geo_so$day.celltype)

# Consider pericyte cluster D21 v D7 --------------------------------------

# Reset cell identities to the combined condition + cluster label
Idents(geo_so) = 'day.celltype'

# run comparison for D21 vs D0, using wilcoxon test
de_cell_pericyte_D21_vs_D7 = FindMarkers(
  object = geo_so,
  slot = 'data', test = 'wilcox',
  ident.1 = 'Day21_Pericyte', ident.2 = 'Day7_Pericyte')

head(de_cell_pericyte_D21_vs_D7)

# Add gene symmbols names and save ----------------------------------------

# Add rownames as a column for output
de_cell_pericyte_D21_vs_D7$gene = rownames(de_cell_pericyte_D21_vs_D7)

# save
write_csv(de_cell_pericyte_D21_vs_D7, 
          file = 'results/tables/de_standard_pericyte_D21_vs_D7.csv')

# summarize diffex results
table(de_cell_pericyte_D21_vs_D7$p_val_adj < 0.05 & 
        abs(de_cell_pericyte_D21_vs_D7$avg_log2FC) > 1.5)

# Create pseudobulk object -------------------------------------------------
pseudo_catch_so = AggregateExpression(geo_so, assays = 'RNA', return.seurat = TRUE, group.by = c('cell_type', 'day', 'replicate'))

# Set up labels to use for comparisons & assign as cell identities
pseudo_catch_so$day.celltype = paste(pseudo_catch_so$day, pseudo_catch_so$cell_type, sep = '_')
Idents(pseudo_catch_so) = 'day.celltype'


# Run pseudobulk comparison between Day 21 and Day 0, using DESeq2 ----------
de_pseudo_pericyte_D21_vs_D7 = FindMarkers(
  object = pseudo_catch_so, 
  ident.1 = 'Day21_Pericyte', ident.2 = 'Day7_Pericyte', 
  test.use = 'DESeq2')

# Take a look at the table
head(de_pseudo_pericyte_D21_vs_D7)

# Add genes and review pseudobulk results ----------------------------------

# Add rownames as a column for output
de_pseudo_pericyte_D21_vs_D7$gene = rownames(de_pseudo_pericyte_D21_vs_D7)

# save results
write_csv(de_pseudo_pericyte_D21_vs_D7,
          file = 'results/tables/de_pseudo_pericyte_D21_vs_D7.csv')

# look at results, using the same thresholds
table(de_pseudo_pericyte_D21_vs_D7$p_val_adj < 0.05 & 
        abs(de_pseudo_pericyte_D21_vs_D7$avg_log2FC) > 1.5)

# Make a volcano plot of pseudobulk diffex results ------------------------
pseudo_pericyte_D21_vs_D7_volcano = ggplot(de_pseudo_pericyte_D21_vs_D7, aes(x = avg_log2FC, y = -log10(p_val))) + geom_point()

ggsave(filename = 'results/figures/volcano_de_pseudo_pericyte_D21_vs_D0.png', 
       plot = pseudo_pericyte_D21_vs_D7_volcano, width = 7, height = 7, units = 'in')

pseudo_pericyte_D21_vs_D7_volcano


# UMAP feature plot of Cd55 gene ------------------------------------------
FeaturePlot(geo_so, features = "Cd55", split.by = "day")
#FeaturePlot(geo_so, features = "Cd55", split.by = "day", label = TRUE)

# Discard all ggplot objects currently in environment ---------------------
# Ok since we saved the plots as we went along
rm(list=names(which(unlist(eapply(.GlobalEnv, is.ggplot))))); 
gc()

# Save Seurat object ------------------------------------------------------
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated_final.rds')

