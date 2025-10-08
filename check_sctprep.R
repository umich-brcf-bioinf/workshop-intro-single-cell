library(here)
library(Seurat)
library(tidyverse)
library(Matrix)
library(BPCells)

# -------------------------------------------------------------------------
# Create seurat object

sample_dirs = list(
  HODay0replicate1  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate1", 
  HODay0replicate2  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate2",
  HODay0replicate3  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate3",
  HODay0replicate4  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay0replicate4",
  HODay7replicate1  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate1",
  HODay7replicate2  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate2",
  HODay7replicate3  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate3",
  HODay7replicate4  = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay7replicate4",
  HODay21replicate1 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate1",
  HODay21replicate2 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate2",
  HODay21replicate3 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate3",
  HODay21replicate4 = "~/ISC_R/inputs/10x_cellranger_filtered_triples/count_run_HODay21replicate4")

# Create the expression matrix from sample dirs
#   Read10X needs a *vector* instead of a *list*, so we use *unlist* to convert
geo_mat = Read10X(data.dir = unlist(sample_dirs))

#geo_mat = open_matrix_dir(dir = '~/ISC_R/bpcells')
geo_so = CreateSeuratObject(counts = geo_mat, min.cells = 1, min.features = 50)
geo_so
str(GetAssayData(geo_so, layer='counts')@x)


summary(GetAssayData(geo_so, layer='counts')[0:4000, 0:4000]@x)

summary(GetAssayData(geo_so, layer='counts')[0:4000, 0:4000]@p)

data = data.frame(when='before_normalization', source = geo, value = GetAssayData(geo_so, layer='counts')[0:4000, 0:4000]@x)

ggplot(data, aes(x = value)) + 
  geom_histogram() +
  scale_y_log10(oob = scales::squish_infinite)

geo_so[['RNA']] = split(geo_so[['RNA']], f = geo_so$orig.ident)
geo_so

data = rbind(data.frame(when='before_normalization', source = 'counts.HODay0replicate1', value = GetAssayData(geo_so, layer='counts.HODay0replicate1')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay0replicate2', value = GetAssayData(geo_so, layer='counts.HODay0replicate2')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay0replicate3', value = GetAssayData(geo_so, layer='counts.HODay0replicate3')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay0replicate4', value = GetAssayData(geo_so, layer='counts.HODay0replicate4')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay7replicate1', value = GetAssayData(geo_so, layer='counts.HODay7replicate1')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay7replicate2', value = GetAssayData(geo_so, layer='counts.HODay7replicate2')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay7replicate3', value = GetAssayData(geo_so, layer='counts.HODay7replicate3')[0:500, 0:500]@x),
             data.frame(when='before_normalization', source = 'counts.HODay7replicate4', value = GetAssayData(geo_so, layer='counts.HODay7replicate4')[0:500, 0:500]@x),
              data.frame(when='before_normalization', source = 'counts.HODay21replicate1', value = GetAssayData(geo_so, layer='counts.HODay21replicate1')[0:500, 0:500]@x),
              data.frame(when='before_normalization', source = 'counts.HODay21replicate2', value = GetAssayData(geo_so, layer='counts.HODay21replicate2')[0:500, 0:500]@x),
              data.frame(when='before_normalization', source = 'counts.HODay21replicate3', value = GetAssayData(geo_so, layer='counts.HODay21replicate3')[0:500, 0:500]@x),
              data.frame(when='before_normalization', source = 'counts.HODay21replicate4', value = GetAssayData(geo_so, layer='counts.HODay21replicate4')[0:500, 0:500]@x))

ggplot(data, 
       aes(x = value, color=source, group=source)) + 
  geom_density(bw=150) +
  xlim(0,500) +
  scale_y_log10()


ggplot(data, 
       aes(x=source, y = value, color=source, group=source)) + 
  geom_boxplot() +
  scale_y_log10()

so = readRDS(here('source/results/rdata/geo_so_sct_normalized.rds'))
so



#Assays(so)
#Layers(so)

#before_counts = GetAssayData(so, layer='counts')
#before_data = GetAssayData(so, layer='data')
before_scale_data = GetAssayData(so, layer='scale.data')
ncol(before_scale_data)

before_scale_data[0:5, 0:5]

# Create a sample dgCMatrix

# hist(as(before_counts, "TsparseMatrix")@x)
# hist(as(before_data, "TsparseMatrix")@x)
# hist(as(before_scale_data, "TsparseMatrix")@x)


data = rbind(data.frame(when='before', source = 'counts', value = as(GetAssayData(so, layer='counts')[0:4000, 0:4000], "TsparseMatrix")@x),
             data.frame(when='before', source = 'data', value = as(GetAssayData(so, layer='data')[0:4000, 0:4000], "TsparseMatrix")@x),
             data.frame(when='before', source = 'scale_data', value = as(GetAssayData(so, layer='scale.data')[0:4000, 0:4000], "TsparseMatrix")@x))



summary(as(GetAssayData(so, layer='counts')[0:4000, 0:4000], "TsparseMatrix")@x)
summary(as(GetAssayData(so, layer='counts')[0:4000, 0:4000], "TsparseMatrix")@x)
summary(as(GetAssayData(so, layer='scale.data')[0:4000, 0:4000], "TsparseMatrix")@x)

ggplot(data, aes(x = value, group = source)) + 
  geom_histogram() +
  scale_y_log10(oob = scales::squish_infinite) +
  facet_grid(. ~ source, scales = "free_x")



so = readRDS(here('source/results/rdata/geo_so_sct_integrated_with_markers.rds'))
so
#Assays(so)
#Layers(so)

#before_counts = GetAssayData(so, layer='counts')
#before_data = GetAssayData(so, layer='data')

# Create a sample dgCMatrix

# hist(as(before_counts, "TsparseMatrix")@x)
# hist(as(before_data, "TsparseMatrix")@x)
# hist(as(before_scale_data, "TsparseMatrix")@x)


data = rbind(data.frame(when='after', source = 'counts', value = as(GetAssayData(so, layer='counts')[0:4000, 0:4000], "TsparseMatrix")@x),
             data.frame(when='after', source = 'data', value = as(GetAssayData(so, layer='data')[0:4000, 0:4000], "TsparseMatrix")@x),
             data.frame(when='after', source = 'scale_data', value = as(GetAssayData(so, layer='scale.data')[0:4000, 0:4000], "TsparseMatrix")@x))


# plot 1 as described in question
ggplot(data, aes(x = value, group = source)) + 
  geom_histogram() +
  scale_y_log10(oob = scales::squish_infinite) +
  facet_grid(. ~ source, scales = "free_x")



geo_so = readRDS('/efs/workshop/isc/workshop_setup/make_rdata_relative/rdata-relative/geo_so_unfiltered.rds')
geo_so

setwd('~/ISC_R')

layer1 = GetAssayData(geo_so, layer = 'counts.HODay0replicate1')
dense_count_matrix = as.matrix(layer1)
dim(dense_count_matrix)
genes_measured_per_cell = colSums(dense_count_matrix !=0)

head(genes_measured_per_cell)
summary(genes_measured_per_cell)

# focus on 3rd quartile cells
excerpt_cells = dense_count_matrix[, which(genes_measured_per_cell > 2750 & genes_measured_per_cell < 2850)]
dim(excerpt_cells)

non_empty_genes = excerpt_cells[rowSums(excerpt_cells) > 0, ]
dim(non_empty_genes)
non_empty_genes[0:5, 0:5]

cells_measured_per_gene = rowSums(non_empty_genes != 0) 
summary(cells_measured_per_gene)

excerpt_genes = non_empty_genes[which(cells_measured_per_gene == 37),]
dim(excerpt_genes)

excerpt_genes

tibble(excerpt_genes, )

long_genes = as_tibble(excerpt_genes, rownames = 'gene') %>% 
  pivot_longer(cols = -gene,
               names_to = 'cell',
               values_to='count')
long_genes

ggplot(long_genes, aes(x = cell, y=count)) +
  geom_boxplot() +
  scale_y_log10() + 
  theme(axis.text.x = element_blank())

geo_so = readRDS('/efs/workshop/isc/workshop_setup/make_rdata_relative/rdata-relative/geo_so_sct_normalized.rds')
geo_so

layer2_counts = GetAssayData(geo_so, layer = 'counts')
layer2
normalized_counts = as.matrix(layer2_counts)[unique(long_genes$gene), unique(long_genes$cell)]
dim(normalized_counts)

long_normalized_counts = as_tibble(normalized_counts, rownames = 'gene') %>% 
  pivot_longer(cols = -gene,
               names_to = 'cell',
               values_to='normalized_count')

long_pre_post_normalization = long_genes %>% 
  left_join(long_normalized_counts, by=c('gene', 'cell'))

ggplot(long_pre_post_normalization, aes(x = cell, y=normalized_count)) +
  geom_boxplot() +
  scale_y_log10() + 
  theme(axis.text.x = element_blank())


layer2_scaledata = GetAssayData(geo_so, layer = 'scale.data')
layer2_scaledata
dim(layer2_scaledata)
normalized_scaledata = as.matrix(layer2_scaledata)[,unique(long_genes$cell)]
dim(normalized_scaledata)
head(normalized_scaledata)

measured_scaledata = colSums(normalized_scaledata != 0)
head(measured_scaledata)


geo_so = readRDS('/efs/workshop/isc/workshop_setup/make_rdata_relative/rdata-relative/geo_so_sct_integrated_with_markers.rds')
gc()
geo_so

layer3_counts = GetAssayData(geo_so, layer = 'counts')
layer3_counts
dim(layer3_counts)


long_synthetic_counts = as_tibble(layer3_counts[,unique(long_genes$cell)], rownames = 'gene') %>% 
  pivot_longer(cols = -gene,
               names_to = 'cell',
               values_to='synthetic_count')


long_123_normalization = long_genes %>% 
  left_join(long_normalized_counts, by=c('gene', 'cell')) %>% 
  left_join(long_synthetic_counts, by=c('gene', 'cell'))

head(long_123_normalization)

ggplot(long_123_normalization, aes(x = cell, y=synthetic_count)) +
  geom_boxplot() +
  scale_y_log10() + 
  theme(axis.text.x = element_blank())

geo_so

length(geo_so[["SCT"]]@var.features)


######


# Install if needed
# install.packages(c("glmGamPoi", "ggplot2"))
library(glmGamPoi)
library(ggplot2)

set.seed(42)

# ---- 1. Simulate sequencing depth for 500 cells ----
n_cells <- 500
depth <- rlnorm(n_cells, meanlog = 9, sdlog = 0.6)  # total UMIs per cell

# ---- 2. Simulate counts for one gene ----
# Expected mean expression increases ~linearly with log(depth)
true_beta0 <- -5.5    # intercept
true_beta1 <- 1.0     # slope vs log(depth)
mu <- exp(true_beta0 + true_beta1 * log(depth))
theta <- 10           # dispersion parameter (larger = less overdispersion)

# Generate observed UMI counts from NB
counts <- rnbinom(n_cells, mu = mu, size = theta)

# ---- 3. Fit a negative binomial regression ----
fit <- glm_gp(
  counts ~ log(depth),
  overdispersion = TRUE,
  size_factors = FALSE
)

summary(fit)

# Extract fitted mean and residuals
#fitted_mu <- fitted(fit)
# If fitted(fit) returns NULL, use this:
beta <- as.numeric(fit$Beta)  # intercept and slope coefficients
eta <- beta[1] + beta[2] * log(depth)
fitted_mu <- exp(eta)

summary(fit)

pearson_resid = (counts - fitted_mu) / sqrt(fitted_mu + fitted_mu^2 / fit$overdispersions)

# ---- 4. Plot observed vs fitted ----
df <- data.frame(
  depth = depth,
  log_depth = log(depth),
  counts = counts,
  fitted_mu = fitted_mu,
  pearson_resid = pearson_resid
)

# (A) Relationship between counts and depth
p1 <- ggplot(df, aes(x = log_depth, y = counts)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = fitted_mu), color = "red", linewidth = 1) +
  labs(
    title = "Negative binomial fit (as used in SCTransform)",
    x = "log(Sequencing depth)",
    y = "Observed counts (GeneX)"
  ) +
  theme_bw()

# (B) Pearson residuals (variance-stabilized)
p2 <- ggplot(df, aes(x = log_depth, y = pearson_resid)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_hline(yintercept = 0, color = "red") +
  labs(
    title = "Pearson residuals (depth-corrected)",
    x = "log(Sequencing depth)",
    y = "Residuals"
  ) +
  theme_bw()

p1
p2

####

# Install if needed
# install.packages(c("Seurat", "ggplot2"))
library(Seurat)
library(ggplot2)

set.seed(42)

# Simulate a small dataset (100 genes x 300 cells)
counts <- matrix(rpois(100 * 300, lambda = rlnorm(100 * 300, 1, 1)), nrow = 100)
rownames(counts) <- paste0("Gene", 1:100)
colnames(counts) <- paste0("Cell", 1:300)

# Create a Seurat object
obj <- CreateSeuratObject(counts = counts)

# Run SCTransform (models depth automatically)
obj <- SCTransform(obj, verbose = FALSE)

# Pick one gene to inspect
gene <- "Gene10"

# Extract all three layers
raw_counts <- GetAssayData(obj, assay = "RNA", slot = "counts")[gene, ]
sct_counts <- GetAssayData(obj, assay = "SCT", layer = "counts")[gene, ]
sct_resid  <- GetAssayData(obj, assay = "SCT", slot = "scale.data")[gene, ]
depth      <- colSums(obj[["RNA"]]$counts)

# Combine into one data frame
df <- data.frame(
  cell = names(raw_counts),
  depth = depth,
  RNA_counts = as.numeric(raw_counts),
  SCT_counts = as.numeric(sct_counts),
  SCT_residuals = as.numeric(sct_resid)
)

head(df)



# Plot 1: raw vs SCT counts vs sequencing depth
p1 <- ggplot(df, aes(x = depth)) +
  geom_point(aes(y = RNA_counts, color = "Raw counts"), alpha = 0.6) +
  geom_point(aes(y = SCT_counts, color = "SCT counts"), alpha = 0.6) +
  scale_x_log10() +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = paste0(gene, ": Raw vs SCT counts"),
    x = "Sequencing depth (total UMIs per cell)",
    y = "Expression (log1p scale)",
    color = "Layer"
  ) +
  theme_bw()
library(tidyverse)

df %>% 
  mutate(qtile = ntile(depth, n = 5)) %>% 
  ggplot(aes(x = RNA_counts, y=SCT_counts, color=qtile)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0(gene, ": Raw vs SCT counts"),
    x = "RNA counts",
    y = "SCT counts") +
  theme_bw()

p3


# Plot 2: residuals vs depth
p2 <- ggplot(df, aes(x = depth, y = SCT_residuals)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_log10() +
  labs(
    title = paste0(gene, ": Pearson residuals (variance-stabilized)"),
    x = "Sequencing depth (total UMIs per cell)",
    y = "SCT residuals"
  ) +
  theme_bw()

p1
p2

#####

# install.packages(c("glmGamPoi", "ggplot2"))  # if needed
library(glmGamPoi)
library(ggplot2)

set.seed(1)

# ---- 1) Simulate cells and an NB relationship with depth ----
n_cells <- 400
depth <- rlnorm(n_cells, meanlog = 9, sdlog = 0.5)      # total UMIs per cell (library size)
beta0 <- -5.2                                            # intercept
beta1 <- 1.0                                             # slope vs log(depth)
mu_true <- exp(beta0 + beta1 * log(depth))               # true mean
theta <- 12                                              # NB size (larger -> less overdispersion)
y <- rnbinom(n_cells, mu = mu_true, size = theta)        # observed RNA counts (UMIs) for GeneX

# ---- 2) Fit a gamma-Poisson / NB regression (as in SCTransform) ----
fit <- glm_gp(y ~ log(depth), overdispersion = TRUE, size_factors = FALSE, on_disk = FALSE)

# Fitted means (if fitted() returns NULL on some versions, compute manually)
fitted_mu <- tryCatch(fitted(fit), error = function(e) NULL)
if (is.null(fitted_mu)) {
  b <- as.numeric(fit$Beta)
  fitted_mu <- exp(b[1] + b[2] * log(depth))
}

# Extract dispersion (NB parameterization used below)
phi <- as.numeric(fit$overdispersions)  # gamma-Poisson overdispersion
# For Pearson residuals we need size = 1/phi; variance = mu + mu^2 * phi

# ---- 3) Compute SCT "adjusted counts" at a reference depth ----
# SCTransform conceptually reports predicted counts as if all cells had the same depth.
# A common choice is the median depth (you could choose mean or any fixed reference).
ref_depth <- median(depth)
b <- as.numeric(fit$Beta)
# Predicted mean at the reference depth for *every* cell:
sct_counts <- exp(b[1] + b[2] * log(ref_depth))

# Note: sct_counts is a single number per gene replicated for all cells.
# To keep shapes consistent, repeat to length n_cells:
sct_counts <- rep(sct_counts, length.out = n_cells)

# ---- 4) Compute Pearson residuals (what SCT stores in scale.data) ----
pearson_resid <- (y - fitted_mu) / sqrt(fitted_mu + (fitted_mu^2) * phi)

# ---- 5) Visualize: raw vs SCT counts; residuals vs depth ----
df <- data.frame(
  depth = depth,
  RNA_counts = y,
  Fitted_mu = fitted_mu,
  SCT_counts = sct_counts,
  Residuals = pearson_resid
)

# (A) Raw counts depend on depth; SCT adjusted counts are constant at ref depth
p1 <- ggplot(df, aes(x = log(depth))) +
  geom_point(aes(y = RNA_counts), alpha = 0.5) +
  geom_line(aes(y = Fitted_mu), color = "red") +
  geom_hline(aes(yintercept = unique(SCT_counts)), linetype = 2) +
  labs(
    title = "Raw counts vs fitted mean and SCT adjusted counts",
    x = "log(Sequencing depth)",
    y = "Counts"
  )

# (B) Pearson residuals are centered ~0 with no depth trend
p2 <- ggplot(df, aes(x = log(depth), y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red") +
  labs(
    title = "Pearson residuals (variance-stabilized)",
    x = "log(Sequencing depth)",
    y = "Residual"
  )

p1
p2


ggplot(df, aes(x = log(depth))) +
  geom_point(aes(y = SCT_counts), alpha = 0.5) +
  geom_line(aes(y = Fitted_mu), color = "red") +
  labs(
    title = "Raw counts vs fitted mean and SCT adjusted counts",
    x = "log(Sequencing depth)",
    y = "SCT Counts"
  )

###
geo_so = readRDS('/efs/workshop/isc/workshop_setup/make_rdata_relative/rdata-relative/geo_so_sct_normalized.rds')

head(geo_so[["SCT"]]@var.features, n=20)

geo_so# Extract all three layers
setwd('~/ISC_R')
gene = 'Col1a1'
raw_counts <- GetAssayData(geo_so, assay = "RNA", layer = "counts.HODay0replicate1")[gene, ]
sct_counts <- GetAssayData(geo_so, assay = "SCT", layer = "counts")[gene, colnames(raw_counts)]
sct_resid  <- GetAssayData(geo_so, assay = "SCT", slot = "scale.data")[gene, colnames(raw_counts)]
depth      <- colSums(GetAssayData(geo_so, assay = "RNA", layer = "counts.HODay0replicate1")@matrix[gene,])

length(depth[depth>0])

#head(as.matrix(raw_counts))
#dim(raw_counts)
# Combine into one data frame
df <- data.frame(
  cell = colnames(raw_counts),
  depth = depth,
  RNA_counts = as.numeric(as.matrix(raw_counts)),
  SCT_counts = as.numeric(sct_counts),
  SCT_residuals = as.numeric(sct_resid)
)

head(df)

df %>% 
  mutate(qtile = ntile(depth, n = 5)) %>% 
  ggplot(aes(x = RNA_counts, y=SCT_counts, color=depth)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0(gene, ": Raw vs SCT counts"),
       x = "RNA counts",
       y = "SCT counts") +
  theme_bw() + 
  coord_fixed(ratio = 1)

# p1 <- ggplot(df, aes(x = log(depth))) +
#   geom_point(aes(y = RNA_counts), alpha = 0.5) +
# #  geom_line(aes(y = Fitted_mu), color = "red") +
# #  geom_hline(aes(yintercept = unique(SCT_counts)), linetype = 2) +
#   labs(
#     title = "Raw counts vs fitted mean and SCT adjusted counts",
#     x = "log(Sequencing depth)",
#     y = "Counts"
#   )
# 
# p1
# 
# p2 <- ggplot(df, aes(x = log(depth))) +
#   geom_point(aes(y = SCT_counts), alpha = 0.5) +
#   #  geom_line(aes(y = Fitted_mu), color = "red") +
#   #  geom_hline(aes(yintercept = unique(SCT_counts)), linetype = 2) +
#   labs(
#     title = "Raw counts vs fitted mean and SCT adjusted counts",
#     x = "log(Sequencing depth)",
#     y = "SCT Counts"
#   )
# 
# p2


geo_so = readRDS('/efs/workshop/isc/workshop_setup/make_rdata_relative/rdata-relative/geo_so_sct_integrated_with_markers.rds')
sct_counts_post <- GetAssayData(geo_so, assay = "SCT", layer = "counts")[gene, colnames(raw_counts)]

df %>%
  mutate(sct_counts_undone = sct_counts_post) %>% 
  mutate(qtile = ntile(depth, n = 5)) %>% 
  ggplot(aes(x = SCT_counts, y=sct_counts_undone, color=depth)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0(gene, ": Raw vs SCT counts"),
       x = "SCT counts",
       y = "SCT counts undone") +
  theme_bw() + 
  coord_fixed(ratio = 1)


df %>%
  mutate(sct_counts_undone = sct_counts_post) %>% 
  mutate(qtile = ntile(depth, n = 5)) %>% 
  ggplot(aes(x = RNA_counts, y=sct_counts_undone, color=depth)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(title = paste0(gene, ": Raw vs SCT counts"),
       x = "RNA counts",
       y = "SCT counts undone") +
  theme_bw() + 
  coord_fixed(ratio = 1)

