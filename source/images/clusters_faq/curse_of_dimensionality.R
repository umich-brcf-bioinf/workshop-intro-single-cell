# ----------------------------------------------------------
# Creates two plots to demonstrate the curse of dimensionality

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# ----------------------------------------------------------
# Simulate clusters with distinctly different distances between vs within
# and interatively increase the dimensionality to show how
# the notion of distance falls apart in higher dimensions.

set.seed(42)

# Parameters
dims <- c(2, 3, 10, 100, 1000, 10000, 20000)
n_per_cluster <- 3
cluster_offset <- 60     # separation between cluster centers (in 1st dimension)
cluster_sd <- 3          # standard deviation of each cluster
n_replicates <- 100      # number of random runs per dimension

# Function to simulate two clusters and compute mean distances
simulate_two_clusters <- function(dim) {
  # Generate two Gaussian clusters
  cluster1 <- matrix(rnorm(n_per_cluster * dim, mean = 0, sd = cluster_sd),
                     n_per_cluster, dim)
  cluster2 <- matrix(rnorm(n_per_cluster * dim, mean = 0, sd = cluster_sd),
                     n_per_cluster, dim)
  
  # Offset the second cluster in the first coordinate
  cluster2[, 1] <- cluster2[, 1] + cluster_offset
  
  # Combine and label clusters
  all_points <- rbind(cluster1, cluster2)
  clusters <- rep(c("Cluster1", "Cluster2"), each = n_per_cluster)

  # Compute all pairwise distances
  pairwise <- expand_grid(i = 1:nrow(all_points), j = 1:nrow(all_points)) %>%
    filter(i < j) %>%
    mutate(
      dist = sqrt(rowSums((all_points[i, ] - all_points[j, ])^2)),
      same_cluster = clusters[i] == clusters[j]
    )
  
  # Return mean distances
  tibble(
    Mean_within = mean(pairwise$dist[pairwise$same_cluster]),
    Mean_between = mean(pairwise$dist[!pairwise$same_cluster])
  )
}

# Run simulation: 100 replicates per dimension
results <- expand_grid(Dimension = dims, Replicate = 1:n_replicates) %>%
  rowwise() %>%
  mutate(vals = list(simulate_two_clusters(Dimension))) %>%
  unnest(cols = c(vals)) %>%
  group_by(Dimension) %>%
  summarise(
    avg_dist_within_cluster = round(mean(Mean_within),0),
#    SD_within = sd(Mean_within),
    avg_dist_between_cluster = round(mean(Mean_between),0),
#    SD_between = sd(Mean_between),
    ratio = round(100 *mean(Mean_within / Mean_between),0)
  )

# Show final summary
print(results)

results %>%
  rename(between_cluster = avg_dist_between_cluster,
         within_cluster = avg_dist_within_cluster) %>% 
  mutate(dimension_f = factor(Dimension)) %>%
  select(dimension_f, between_cluster, within_cluster) %>% 
  pivot_longer(cols=-dimension_f, names_to='metric', values_to='dist') %>% 
  ggplot(aes(x = dimension_f, y = dist, group = metric, color = metric)) +
  geom_point(size=3, alpha=0.6) +
  geom_line() +
  theme_minimal(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.title = element_blank(), 
        legend.position = 'top') +
  labs(y='distance', x='dimensions') +
  scale_y_log10()

ggsave(filename=here('source/images/clusters_faq/curse_of_dimensionality.png'), 
       height = 3, width = 4, units = 'in')

# ----------------------------------------------------------
# Simulate a 2D cluster with distinctly different distances between vs within

#----------------
set.seed(42)

# Create two clusters of 3 points each in 2D
cluster1 <- data.frame(
  x = rnorm(3, mean = 20, sd = 3),
  y = rnorm(3, mean = 20, sd = 3),
  Cluster = "Cluster 1"
)

cluster2 <- data.frame(
  x = rnorm(3, mean = 80, sd = 3),
  y = rnorm(3, mean = 80, sd = 3),
  Cluster = "Cluster 2"
)

# Combine into one dataset
data <- rbind(cluster1, cluster2)

data$shape = ifelse(data$Cluster=='Cluster 1', 17, 15)
data
# Plot with ggplot2
p = ggplot(data, aes(x = x, y = y, color = Cluster, shape=shape)) +
  geom_point(size = 3) +
  scale_shape_identity() +
  coord_equal() +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text=element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 1))
p

ggsave(filename=here('source/images/clusters_faq/curse_of_dimensionality_cluster2d.png'), 
  height = 1.5, width = 1.5, units = 'in')
