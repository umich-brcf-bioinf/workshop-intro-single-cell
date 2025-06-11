library(Seurat)
library(tidyverse)

geo_so = readRDS('~/ISC_R/inputs/prepared_data/rdata/geo_so_sct_clustered.rds')

# -------------------------------------------------------------------------
resolutions = c(0.2, 0.4, 0.8, 1.2, 1.6)
included_pcs = c(5, 10 , 25, 50)

# get_data = function (object, group.by, reduction, dims = c(1,2)) {
#   cells <- Cells(x = object, assay = DefaultAssay(object = object[[reduction]]))
#   dims <- paste0(Key(object = object[[reduction]]), dims)
#   orig.groups <- group.by
#   group.by <- group.by %||% 'ident'
#   data <- FetchData(
#     object = object,
#     vars = c(dims, group.by),
#     cells = cells,
#     clean = 'project')
#   return (data)
# }

pc_res_cluster = tibble(num_pcs = numeric(), resolution= numeric(), cluster_count=numeric())

for (pc in included_pcs) { 
  for(res in resolutions) {
    name_suffix = sprintf('pcs=%s_res=%s', pc, res)
    umap_file = sprintf('results/figures/cluster_matrix_%s.png', name_suffix)
    cluster_column = sprintf('SCT_snn_res.%s', res)
    
    message(name_suffix)
    
    geo_so = suppressMessages(FindNeighbors(geo_so, dims = 1:pc, reduction = 'integrated.sct.rpca'))
    geo_so = suppressMessages(FindClusters(geo_so, resolution = res))
    cluster_count = length(unique(geo_so@meta.data[[cluster_column]]))
    
    pc_res_cluster = pc_res_cluster %>%
      add_row(num_pcs = pc, resolution = res, cluster_count=cluster_count)
    
    # DimPlot(geo_so,
    #         group.by = cluster_column,
    #         label = TRUE,
    #         reduction = 'umap.integrated.sct.rpca') +
    #   NoLegend() +
    #   labs(title=paste(name_suffix, ':', cluster_count, 'clusters')) +
    #   theme(axis.title=element_blank(),
    #         axis.text=element_blank(),
    #         axis.ticks=element_blank(),
    #         axis.line = element_blank())
    # 
    # ggsave(filename = umap_file,
    #        width = 8, height = 7, units = 'in')
    
    # data = get_data(geo_so, 
    #                 group.by = cluster_column,
    #                 reduction = 'umap.integrated.sct.rpca')
    # write_csv(data, file = sprintf('results/tables/%s.csv', name_suffix))
  }
}

write_csv(pc_res_cluster, '~/pc_res_cluster.csv')

pc_res_cluster$num_pcs = factor(pc_res_cluster$num_pcs)

pc_res_cluster %>% 
  ggplot(aes(x=resolution, y=cluster_count, group=num_pcs, color=num_pcs)) +
    geom_line() +
    theme_minimal()
  