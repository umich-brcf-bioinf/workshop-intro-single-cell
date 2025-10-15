# =========================================================================
# Day 3 set up
# =========================================================================

# -------------------------------------------------------------------------
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
options(future.globals.maxSize = 1e9)


# -------------------------------------------------------------------------
# Set number of principal components
pcs = 10

# -------------------------------------------------------------------------
# Load seurat object
geo_so = readRDS('~/ISC_R/inputs/prepared_data/rdata/geo_so_sct_clustered.rds')
geo_so


# =========================================================================
# Marker identification and visualization
# =========================================================================

# Prep for cluster comparisons to find empirical markers
geo_so = SetIdent(geo_so, value = 'integrated.sct.rpca.clusters')
?PrepSCTFindMarkers
geo_so = PrepSCTFindMarkers(geo_so)

# Run comparisons for each cluster to generate markers
geo_markers = FindAllMarkers(geo_so, only.pos = TRUE, min.pct = 0.05)

# Write out full cluster marker results to file
write_csv(geo_markers, file = 'results/tables/marker_genes_PC10-0.4res.csv')

# Take a look at the first few rows of the result
head(geo_markers)

# -------------------------------------------------------------------------
# Identify marker genes for each cluster
# Create table of top 5 markers per cluster (using default ranking)
top_5 = geo_markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% slice_head(n = 5)

# Look at results
head(top_5, n = 10)

# -------------------------------------------------------------------------
# Visualize top marker genes as dot plot
top_5_sct_dot_plot = DotPlot(geo_so, features = unique(top_5$gene)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = 'Top 5 Cluster Genes by FDR and avg_log2FC') + coord_flip()
top_5_sct_dot_plot

# -------------------------------------------------------------------------
# Save dot plot of top marker genes
ggsave(filename = 'results/figures/markers_top_5_sct_dot_plot.png', 
       plot = top_5_sct_dot_plot, 
       width = 8, height = 18, units = 'in') 

# -------------------------------------------------------------------------
# Check mitochondrial gene expression
percent_mito_plot = FeaturePlot(geo_so, features='percent.mt')
percent_mito_plot


# ---------------------------
# ADDED - check an example marker gene
FeaturePlot(geo_so, features='Acod1')


# save to file
ggsave(filename = 'results/figures/percent_umap_mito_plot.png', 
       plot = percent_mito_plot, 
       width = 6, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Remove the plot variables from the environment to avoid excessive memory usage
plots = c("top_5_sct_dot_plot", 
          "top_5_rna_dot_plot", 
          "percent_mito_plot")

# Only remove plots that actually exist in the environment
rm(list=Filter(exists, plots))
gc()

# -------------------------------------------------------------------------
# Save Seurat object and gene marker data
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated_with_markers.rds')
saveRDS(geo_markers, file = 'results/rdata/geo_markers.rds')

# =========================================================================
# Cell Type Annotation
# =========================================================================

# Load scCATCH
library(scCATCH)

# check that cell identities are set to expected resolution 
all(Idents(geo_so) == geo_so$integrated.sct.rpca.clusters)

# -------------------------------------------------------------------------
# Annotate clusters using scCATCH

# create scCATCH object, using count data
geo_catch = createscCATCH(data = geo_so@assays$SCT@counts, cluster = as.character(Idents(geo_so)))

# add marker genes to use for predictions
geo_catch@markergene = geo_markers

# specify tissues/cell-types from the scCATCH reference
tissues_of_interest = c('Blood', 'Peripheral Blood', 'Muscle', 'Skeletal muscle', 'Epidermis', 'Skin')
species_tissue_mask = cellmatch$species == 'Mouse' & cellmatch$tissue %in% tissues_of_interest
geo_catch@marker = cellmatch[species_tissue_mask, ]

# run scCATCH to generate predictions
?findcelltype
geo_catch = findcelltype(geo_catch)
head(geo_catch) # whole object
head(geo_catch@celltype) # cell type prediction slot

# look at the predictions
sc_pred_tbl = geo_catch@celltype %>% 
  select(cluster, cell_type, celltype_score) %>%
  rename(cell_type_prediction = cell_type)
sc_pred_tbl

# -------------------------------------------------------------------------
# Write scCATCH predictions to file
write_csv(sc_pred_tbl, file = 'results/tables/scCATCH_cluster_predictions.csv')

# -------------------------------------------------------------------------
# Create lists of immune cells and associated gene markers
immune_markers = list()
immune_markers[['Inflam. Macrophage']] = c('Cd14', 'Cxcl2') # Cd14 a- monocyte/macrophage cells
immune_markers[['Platelet']] = c('Pf4')
immune_markers[['Mast cells']] = c('Gata2', 'Kit')
immune_markers[['NK cells']] = c('Nkg7', 'Klrd1')
immune_markers[['B-cell']] = c( 'Ly6d', 'Cd19', 'Cd79b', 'Ms4a1')
immune_markers[['T-cell']] = c( 'Cd3d','Cd3e','Cd3g') # also Thy1


# -------------------------------------------------------------------------
# Plot other immune to assist with cluster identification
immune_markers_plot = DotPlot(geo_so, features = immune_markers, assay = 'SCT')  +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 45, vjust = 0.5))
immune_markers_plot

# save to file
ggsave(filename = 'results/figures/immune_markers_sct_dot_plot.png', 
       plot = immune_markers_plot, width = 10, height = 5, units = 'in')


# -------------------------------------------------------------------------
# Create lists of other cells and associated gene markers
other_markers = list()
other_markers[['Pericyte']] = c('Acan','Sox9')
other_markers[['SMC']] = c('Acta2', 'Myh11') # SMC = mesenchymal smooth-muscle cell/mesenchymal lineage
other_markers[['Keratinocytes']] = c('Thy1', 'Dlk1') # fibro progenitors aso=Thy1
other_markers[['Myofibroblasts']] = c('Tmem100', 'Cd34', 'Ly6c1') # hematopoetic stem/activated fibroblast=Cd34
other_markers[['Fibroblast']] = c('Dpt', 'Fn1', 'Col3a1')  # activated fib = Fn1
other_markers[['Endothelial']] = c('Pecam1', 'Cd38') # from wound healing; Pecam1 also exp in endothelial
other_markers[['HSC']] = c('Ltb', 'Cd74') # less well defined/conflicting definitions
other_markers[['Erythroid']] = c('Hba-a1')

# -------------------------------------------------------------------------
# Plot known cell-type markers
other_markers_dot_plot = DotPlot(geo_so, features = other_markers, assay = 'SCT') +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 45, vjust = 0.5))
other_markers_dot_plot

# save to file
ggsave(filename = 'results/figures/other_markers_sct_dot_plot.png', 
       plot = other_markers_dot_plot, width = 12, height = 5, units = 'in')

# NOT IN SHOWNOTES - how to add labels for long version of DotPlot()
immune_markers_plot + coord_flip() +
  facet_grid(rows = vars(feature.groups),
             scales = 'free', space = 'free', switch ='y') +
  theme(strip.placement = 'outside',
        strip.text.y.left = element_text(angle = 0))

# -------------------------------------------------------------------------
# Load the revised annotations from prepared file
celltype_annos = read_csv('inputs/prepared_data/revised_cluster_annotations.csv') %>%
  mutate(cluster=factor(cluster))
head(celltype_annos)
View(celltype_annos)

# -------------------------------------------------------------------------
# Merge cell types in but as a new table to slide into @meta.data
copy_metadata = geo_so@meta.data
new_metadata = copy_metadata %>% 
  left_join(celltype_annos, by = c('integrated.sct.rpca.clusters' = 'cluster'))
#  We are implicitly relying on the same row order!
rownames(new_metadata) = rownames(geo_so@meta.data)
head(new_metadata)

# Replace the meta.data
geo_so@meta.data = new_metadata 

head(geo_so@meta.data)

# -------------------------------------------------------------------------
# Make a labeled UMAP plot of clusters
annos_umap_plot = 
  DimPlot(geo_so, group.by = 'cell_type', label = TRUE, reduction = 'umap.integrated.sct.rpca')
annos_umap_plot

ggsave(filename = 'results/figures/umap_integrated_annotated.png', 
       plot =  annos_umap_plot,
       width = 8, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Repeat, splitting by day
annos_umap_condition_plot = 
  DimPlot(geo_so,
          group.by = 'cell_type',
          split.by = 'time',
          label = FALSE,
          reduction = 'umap.integrated.sct.rpca')
annos_umap_condition_plot

ggsave(filename = 'results/figures/umap_integrated_annotated_byCondition.png', 
       plot = annos_umap_condition_plot,
       width = 14, height = 5, units = 'in')

# -------------------------------------------------------------------------
# Discard all ggplot objects currently in environment
# (Ok since we saved the plots as we went along)
rm(list=names(which(unlist(eapply(.GlobalEnv, is.ggplot))))); 
gc()


# -------------------------------------------------------------------------
# Save Seurat object and annotations
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated_with_annos.rds')
saveRDS(geo_catch, file = 'results/rdata/geo_catch.rds')

# =========================================================================
# Differential Expression Analysis
# =========================================================================

# Create combined label of day + celltype to use for DE contrasts
geo_so$day.celltype = paste(geo_so$time, geo_so$cell_type, sep = '_')
# check labels
unique(geo_so$day.celltype)

# Reset cell identities to the combined condition + celltype label
head(Idents(geo_so))
Idents(geo_so) = 'day.celltype'
head(Idents(geo_so))

# -------------------------------------------------------------------------
# Consider pericyte cluster D21 v D7 & run DE comparison using wilcoxon test
de_cell_pericyte_D21_vs_D7 = FindMarkers(
  object = geo_so,
  slot = 'data', test = 'wilcox',
  ident.1 = 'Day21_Pericyte', ident.2 = 'Day7_Pericyte')

head(de_cell_pericyte_D21_vs_D7)

# -------------------------------------------------------------------------
# Add gene symbols names and save

# Add rownames as a column for output
de_cell_pericyte_D21_vs_D7$gene = rownames(de_cell_pericyte_D21_vs_D7)

# save to file
write_csv(de_cell_pericyte_D21_vs_D7, 
          file = 'results/tables/de_standard_pericyte_D21_vs_D7.csv')

# summarize diffex results
table(de_cell_pericyte_D21_vs_D7$p_val_adj < 0.05 & 
        abs(de_cell_pericyte_D21_vs_D7$avg_log2FC) > 1.5)

# -------------------------------------------------------------------------
# Compare pericyte cluster D7 v D0
de_cell_pericyte_D7_vs_D0 = FindMarkers(
  object = geo_so,
  slot = 'data', test = 'wilcox',
  ident.1 = 'Day7_Pericyte', ident.2 = 'Day0_Pericyte')

head(de_cell_pericyte_D7_vs_D0)

# -------------------------------------------------------------------------
# Add rownames for D7 v D0 results
de_cell_pericyte_D7_vs_D0$gene = rownames(de_cell_pericyte_D7_vs_D0)

# summarize results
table(de_cell_pericyte_D7_vs_D0$p_val_adj < 0.05 & 
        abs(de_cell_pericyte_D7_vs_D0$avg_log2FC) > 1.5)


# -------------------------------------------------------------------------
# Create pseudobulk object
pseudo_catch_so = 
  AggregateExpression(geo_so, 
                      assays = 'RNA',
                      return.seurat = TRUE,
                      group.by = c('cell_type', 'time', 'replicate'))

# Set up labels to use for comparisons & assign as cell identities
pseudo_catch_so$day.celltype = paste(pseudo_catch_so$time, pseudo_catch_so$cell_type, sep = '_')
Idents(pseudo_catch_so) = 'day.celltype'

# -------------------------------------------------------------------------
# Run pseudobulk comparison between Day 21 and Day 0, using DESeq2
de_pseudo_pericyte_D21_vs_D7 = FindMarkers(
  object = pseudo_catch_so, 
  ident.1 = 'Day21_Pericyte', ident.2 = 'Day7_Pericyte', 
  test.use = 'DESeq2')

# Take a look at the table
head(de_pseudo_pericyte_D21_vs_D7)

# -------------------------------------------------------------------------
# add genes rownames as a column for output
de_pseudo_pericyte_D21_vs_D7$gene = rownames(de_pseudo_pericyte_D21_vs_D7)

# save results
write_csv(de_pseudo_pericyte_D21_vs_D7,
          file = 'results/tables/de_pseudo_pericyte_D21_vs_D7.csv')

# review pseudobulk results, using the same thresholds
table(de_pseudo_pericyte_D21_vs_D7$p_val_adj < 0.05 & 
        abs(de_pseudo_pericyte_D21_vs_D7$avg_log2FC) > 1.5)

# DEMO - volcano plot
# -------------------------------------------------------------------------
# Make a volcano plot of pseudobulk diffex results
pseudo_pericyte_D21_vs_D7_volcano = 
  ggplot(de_pseudo_pericyte_D21_vs_D7, aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
  geom_point()
pseudo_pericyte_D21_vs_D7_volcano

ggsave(filename = 'results/figures/volcano_de_pseudo_pericyte_D21_vs_D0.png', 
       plot = pseudo_pericyte_D21_vs_D7_volcano,
       width = 7, height = 7, units = 'in')

# -------------------------------------------------------------------------
# UMAP feature plot of Cd55 gene
FeaturePlot(geo_so, features = "Cd55", split.by = "time")


# -------------------------------------------------------------------------
# Discard all ggplot objects currently in environment
# (Ok since we saved the plots as we went along.)
rm(list=names(which(unlist(eapply(.GlobalEnv, is.ggplot))))); 
gc()

# -------------------------------------------------------------------------
# Save Seurat object
saveRDS(geo_so, file = 'results/rdata/geo_so_sct_integrated_final.rds')
