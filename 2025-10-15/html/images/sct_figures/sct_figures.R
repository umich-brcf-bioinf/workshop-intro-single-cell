library(here)
library(Seurat)
library(BPCells)
library(tidyverse)
library(MASS)

prepared_data_dir = '/efs/workshop/isc/workshop_setup/shared_sample_data/ISC_R/inputs/prepared_data/rdata/'
geo_so = readRDS(here(prepared_data_dir, 'geo_so_sct_normalized.rds'))


median_umi = geo_so[["SCT"]]@SCTModel.list$HODay0replicate1@median_umi
minimum_median_umi = geo_so[["SCT"]]@SCTModel.list$HODay0replicate3@median_umi



# original_gene_umi_scatter_plots -------------------------------------

gene = "C1qb"

# Extract all three layers
rna_counts <- GetAssayData(geo_so, assay = "RNA", layer = "counts.HODay0replicate1")[gene, ]
umi_depth <- colSums(GetAssayData(geo_so, assay = "RNA", layer = "counts.HODay0replicate1")@matrix[,colnames(rna_counts)])
sct_counts <- GetAssayData(geo_so, assay = "SCT", layer = "counts")[gene, colnames(rna_counts)]
sct_data  <- GetAssayData(geo_so, assay = "SCT", layer = "data")[gene, colnames(rna_counts)]
sct_resid  <- GetAssayData(geo_so, assay = "SCT", layer = "scale.data")[gene, colnames(rna_counts)]

#length(umi_depth[umi_depth>0])

#head(as.matrix(raw_counts))
#dim(raw_counts)
# Combine into one data frame
tmp_df <- data.frame(
  cell = colnames(rna_counts),
  umi_count = umi_depth,
  rna_count = as.numeric(as.matrix(rna_counts)),
  sct_count = as.numeric(sct_counts),
  sct_data= as.numeric(sct_data),
  sct_residual = as.numeric(sct_resid)
)

#head(tmp_df)

rna_scatterplot = tmp_df %>% 
  filter(rna_count > 1) %>% 
  ggplot(aes(x = umi_count, y=rna_count)) +
  geom_point(alpha = 0.6, color='blue', size=5) +
  #  geom_smooth(method = "glm", se = FALSE, color = "lightblue") +
  geom_smooth(method = "glm.nb", se = FALSE, color = "lightblue", linewidth=3) +
  geom_vline(xintercept = median_umi, color='red', alpha=0.6, linewidth=2) +
    labs(title = paste0("RNA counts | gene = ", gene, ""),
       x = "log10(cell UMI depth)",
       y = "log10(RNA count)") +
  scale_x_log10(breaks=c(1000, 3000, median_umi, 10000)) +
  scale_y_log10() +
  theme_bw(base_size = 24)

rna_scatterplot
ggsave(filename = here("source/images/sct_figures", "rna_scatterplot.png"),
       plot = rna_scatterplot,
       width = 8, height = 5.5, units = 'in')

sct_counts_scatterplot = tmp_df %>% 
  filter(rna_count > 1) %>% 
  ggplot(aes(x = umi_count, y=sct_count)) +
  geom_point(alpha = 0.6, color='darkgreen', size=5) +
  #geom_smooth(method = "glm", se = FALSE, color = "lightgreen") +
  #geom_smooth(method = "glm.nb", se = FALSE, color = "lightgreen", linewidth=3) +
  #  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0("SCT counts |  gene = ", gene),
       x = "log10(cell UMI depth)",
       y = "SCT count") +
  scale_x_log10() +
  ylim(0,100) +
  theme_bw(base_size = 24)

sct_counts_scatterplot
ggsave(filename = here("source/images/sct_figures", "sct_counts_scatterplot.png"),
       plot = sct_counts_scatterplot,
       width = 8, height = 5.5, units = 'in')

sct_data_scatterplot = tmp_df %>% 
  filter(rna_count > 0) %>% 
  ggplot(aes(x = umi_count, y=sct_data)) +
  geom_point(alpha = 0.6, color='darkgreen', size=5) +
  #geom_smooth(method = "glm", se = FALSE, color = "lightgreen") +
  geom_smooth(method = "glm.nb", se = FALSE, color = "lightgreen", linewidth=3) +
  #  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0("SCT data |  gene = ", gene),
       x = "log10(cell UMI depth)",
       y = "SCT data") +
  scale_x_log10() +
#  ylim(0,100) +
  theme_bw(base_size = 24)

sct_data_scatterplot
ggsave(filename = here("source/images/sct_figures", "sct_data_scatterplot.png"),
       plot = sct_data_scatterplot,
       width = 8, height = 5.5, units = 'in')



sct_scaled_data_scatterplot = tmp_df %>% 
  filter(rna_count > 1) %>% 
  ggplot(aes(x = umi_count, y=sct_residual)) +
  geom_point(alpha = 0.6, color='darkgreen', size=5) +
  #geom_smooth(method = "glm", se = FALSE, color = "lightgreen") +
  # geom_smooth(method = "glm.nb", se = FALSE, color = "lightblue") +
  #  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0("SCT scaled.data | gene = ", gene),
       x = "log10(cell UMI count)",
       y = "scaled.data") +
  scale_x_log10() +
  theme_bw(base_size = 24)

sct_scaled_data_scatterplot
ggsave(filename = here("source/images/sct_figures", "sct_scaled_data_scatterplot.png"),
       plot = sct_scaled_data_scatterplot,
       width = 8, height = 5.5, units = 'in')


# Original boxplots -------------------------------------------------------

variable_genes = geo_so[["SCT"]]@var.features

tmp_df = as.matrix(GetAssayData(geo_so, assay = "SCT", layer = "counts"))[variable_genes,]

df_cell_depth = tibble(
  cell_depth = colSums(tmp_df, na.rm = TRUE),
  cell = colnames(tmp_df),
) %>% 
  mutate(sample = gsub('_.*', '', cell),
         row_id = row_number())

rm(tmp_df)
gc()

original_sct_depth_boxplots = df_cell_depth %>% 
  ggplot(aes(x = sample, y=cell_depth)) +
  geom_boxplot(alpha = 0.6) +
  scale_y_log10() + ylab("log10(cell UMI depth)") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

original_sct_depth_boxplots

ggsave(filename = here("source/images/sct_figures", "original_sct_depth_boxtplots.png"),
       plot = original_sct_depth_boxplots,
       width = 8, height = 5.5, units = 'in')

geo_so2 = readRDS(here(prepared_data_dir, 'geo_so_sct_integrated_with_markers.rds'))
geo_so2 

tmp_df = as.matrix(GetAssayData(geo_so2, assay = "SCT", layer = "counts"))[variable_genes,]

df_cell_depth = tibble(
  cell_depth = colSums(tmp_df, na.rm = TRUE),
  cell = colnames(tmp_df),
) %>% 
  mutate(sample = gsub('_.*', '', cell),
         row_id = row_number())

rm(tmp_df)
gc()

corrected_sct_depth_boxplots = df_cell_depth %>% 
  ggplot(aes(x = sample, y=cell_depth)) +
  geom_boxplot(alpha = 0.6) +
  scale_y_log10() + ylab("log10(cell UMI depth)") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

corrected_sct_depth_boxplots

ggsave(filename = here("source/images/sct_figures", "corrected_sct_depth_boxtplots.png"),
       plot = corrected_sct_depth_boxplots,
       width = 8, height = 5.5, units = 'in')
