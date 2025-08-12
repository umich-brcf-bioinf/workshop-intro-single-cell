# -------------------------------------------------------------------------
# matrix.R
#
# Generates and saves UMAPs and csvs for a matrix of PCs and resolutions
# In addition to the cross-product of resolutions x PCs also generates
#   a "blank" set of blank UMAPs without clusters.
# Takes "a good while" to run
# Overwrites any existing files; does not delete any existing files.
#
# cgates 7/3/2025

here::i_am('source/images/clusters_faq/matrix.R')
library(Seurat)
library(tidyverse)

output_dir = 'source/images/clusters_faq/'
geo_so = readRDS('~/ISC_R/inputs/prepared_data/rdata/geo_so_sct_integrated.rds')

resolutions = c(0.2, 0.4, 0.8, 1.2)
included_pcs = c(5, 10 , 16, 30)

# -------------------------------------------------------------------------

# Initialize a table of stats with the "blank", res=0 results
pc_res_cluster = tibble(expand.grid(included_pcs = included_pcs,
                                    resolution = 0,
                                    cluster_count=0)) %>% 
  mutate(umap_filename=sprintf("umap/res=0_pcs=%02d.png", included_pcs))


# -------------------------------------------------------------------------
# Extracts the cluster data (UMAP1, UMAP2) from a Seurat Object reduction
get_data = function (object, group.by, reduction, dims = c(1,2)) {
  cells <- Cells(x = object, assay = DefaultAssay(object = object[[reduction]]))
  dims <- paste0(Key(object = object[[reduction]]), dims)
  orig.groups <- group.by
  group.by <- group.by %||% 'ident'
  data <- FetchData(
      object = object,
      vars = c(dims, group.by),
      cells = cells,
      clean = 'project') %>% 
    rownames_to_column('cell') %>%
    rename(cluster = group.by)
  return (data)
}



# -------------------------------------------------------------------------
# Run through the combination of PCs and res and generate the UMAPs

for (pc in included_pcs) { 
  for(res in resolutions) {
    name_suffix = sprintf('res=%s_pcs=%02d', res, pc)
    umap_file = sprintf('umaps/%s.png', name_suffix)
    cluster_column = sprintf('SCT_snn_res.%s', res)
    
    message(name_suffix)
    
    geo_so = suppressMessages(FindNeighbors(geo_so, dims = 1:pc, reduction = 'integrated.sct.rpca'))
    geo_so = suppressMessages(FindClusters(geo_so, resolution = res))
    cluster_count = length(unique(geo_so@meta.data[[cluster_column]]))
    
    pc_res_cluster = pc_res_cluster %>%
      add_row(included_pcs = pc, resolution = res, cluster_count=cluster_count, umap_filename=umap_file)
    
    geo_so = RunUMAP(geo_so,
                     dims = 1:pc,
                     reduction = 'integrated.sct.rpca',
                     reduction.name = 'umap.integrated.sct.rpca')
    clustered_umap = DimPlot(geo_so,
            group.by = cluster_column,
            label = TRUE,
            label.size=10,
            reduction = 'umap.integrated.sct.rpca') +
      NoLegend() +
      labs(title=paste(cluster_count, 'clusters')) +
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.line = element_blank())
    ggsave(plot = clustered_umap,
           filename = here::here(output_dir, umap_file),
           width = 8, height = 7, units = 'in', dpi='screen')

    data = get_data(geo_so,
                    group.by = cluster_column,
                    reduction = 'umap.integrated.sct.rpca') %>%
              mutate(res_pc=name_suffix)
    write_csv(data, file = here::here(output_dir, paste0(name_suffix,'.csv.tmp')))

    # Generate a single "blank" umap for each PC
    if (res == min(resolutions)) {
        res = 0
        name_suffix = sprintf('res=%s_pcs=%s', res, pc)
        umap_file = sprintf('umaps/%s.png', name_suffix)

        blank_umap = DimPlot(geo_so,
                         reduction = 'umap.integrated.sct.rpca',
                         cols=rep("#444444", length(unique(Idents(geo_so))))) +
        NoLegend() +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              axis.line = element_blank())

        ggsave(plot = blank_umap,
               filename = here::here(output_dir, umap_file),
               width = 8, height = 7, units = 'in', dpi='screen')
    }
  }
}

files = list.files(path=here::here(output_dir), pattern = "*.csv.tmp", full.names = TRUE)
data = read_csv(files) %>%
  pivot_wider(
    id_cols = "cell",
    names_from = "res_pc",
    values_from = c('umapintegratedsctrpca_1', 'umapintegratedsctrpca_2', 'cluster'),
    names_glue = "{res_pc}|{.value}",
    names_vary='slowest'
  )
write_csv(data, file = here::here(output_dir,'clusters.csv.gz'))
unlink(files)

write_csv(pc_res_cluster, here::here(output_dir, 'pc_res_cluster.csv'))

# -------------------------------------------------------------------------
# generate the line plot of clusters x resolution across PCs

# I know we just wrote this file, but writing it and reading it in separate
# steps makes it easier to tweak this with a new coordinates or revise the 
# aesthetics without having to rerun the whole cross-product or res x PCs.
pc_res_cluster = read_csv(here::here(output_dir, 'pc_res_cluster.csv'))

pc_res_cluster$included_pcs_f = factor(pc_res_cluster$included_pcs)
pc_res_cluster$resolution_f = factor(pc_res_cluster$resolution)


pc_res_cluster_plt = pc_res_cluster %>% 
  ggplot(aes(x= included_pcs_f, y=cluster_count, group=resolution_f, color=resolution_f)) +
  geom_line() +
  geom_point(size=3) +
  theme_minimal(base_size=15) +
  labs(x='Number PCs', y='Number of clusters', color = "Resolution") + 
  guides(color = guide_legend(reverse=TRUE))


pc_res_cluster_plt

ggsave(filename = here::here(output_dir, 'pc_res_cluster.png'),
       plot = pc_res_cluster_plt,
       width = 12, height = 3, units = 'in', dpi='screen')
