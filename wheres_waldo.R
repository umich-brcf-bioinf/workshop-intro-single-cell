# Exact theoretical distance distribution on a 100 x 100 x 100 grid
# (two uniformly random DISTINCT points), then binned for a smooth plot.

library(ggplot2)

n <- 100
d  <- 0:(n - 1)
w1 <- ifelse(d == 0, n, 2 * (n - d))  # ordered-pair weights per offset on one axis
d2 <- d^2

# Build all (dx, dy) squared norms and weights (vectors of length 100*100 = 1e4)
S2_vec <- as.vector(outer(d2, d2, `+`))
W2_vec <- as.vector(outer(w1, w1, `*`))

# Extend to z by adding each dz^2 and multiplying by w1[dz]
# Resulting vectors are length 1e6 (100*100*100)
s_flat <- S2_vec + rep(d2, each = length(S2_vec))
w_flat <- W2_vec * rep(w1, each = length(W2_vec))

# Exclude the zero-offset (dx=dy=dz=0) so pairs are distinct
keep <- s_flat != 0L
s_flat <- s_flat[keep]
w_flat <- w_flat[keep]

# Total mass = number of ordered distinct pairs
N <- n^3
total_mass <- N * (N - 1)

# Aggregate weights by exact squared distance s = dx^2 + dy^2 + dz^2
# (tapply returns a named vector; fill into a complete 0..max_s array)
agg <- tapply(w_flat, s_flat, sum)
max_s <- max(s_flat)
mass_by_s <- numeric(max_s + 1)
mass_by_s[as.integer(names(agg)) + 1] <- as.numeric(agg)

# Convert to probabilities; map s -> r = sqrt(s)
pmf <- mass_by_s / total_mass
support <- which(pmf > 0) - 1L
r <- sqrt(support)
p <- pmf[support + 1L]

# Summary stats
mean_r <- sum(r * p)
sd_r   <- sqrt(sum(p * (r - mean_r)^2))

# Bin radii into small bins to produce a smooth-looking curve
breaks <- seq(0, sqrt(3) * (n - 1), length.out = 50)
bin_id <- cut(r, breaks, include.lowest = TRUE)
bin_prob <- tapply(p, bin_id, sum)
centers <- 0.5 * (head(breaks, -1) + tail(breaks, -1))
widths  <- diff(breaks)
density <- as.numeric(bin_prob) / widths
density[is.na(density)] <- 0

df <- data.frame(distance = centers, density = density)

# Plot
ggplot(df, aes(distance, density)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = mean_r, linetype = "dashed") +
  labs(
    title = "Smoothed theoretical distribution of Euclidean distances",
    subtitle = sprintf("Two random distinct points on a 100×100×100 grid | Mean ≈ %.3f, SD ≈ %.3f", mean_r, sd_r),
    x = "Distance",
    y = "Smoothed probability density"
  ) +
  theme_minimal(base_size = 13)
