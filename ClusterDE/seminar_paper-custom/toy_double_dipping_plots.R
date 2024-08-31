# LLM code
library(mvtnorm)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# Set random seed for reproducibility
set.seed(123)

# Generate multivariate normal data
n <- 1000  # number of points
mu <- c(0, 0)  # mean vector
sigma <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)  # covariance matrix
df <- rmvnorm(n, mean = mu, sigma = sigma) |>
  as.data.frame(make.names = FALSE)
names(df) = c("X1", "X2")

p_mvn_gray <- ggplot(df, aes(x = X1, y = X2)) +
  geom_point(size = 1, shape = 1) +
  labs(title = "Before clustering")

kmeans_result <- kmeans(df, centers = 2)
df$cluster <- as.factor(kmeans_result$cluster)

p_mvn_clustered <- ggplot(df, aes(x = X1, y = X2, color = cluster)) +
  geom_point(size = 1, shape = 1) +
  scale_color_brewer(palette = "Accent") +
  labs(title = "After clustering")

ggarrange(p_mvn_gray, p_mvn_clustered) |>
  ggexport(filename = "../seminar_paper-custom/figures/fig-gray_and_clustered-S15.png",
           width = 800, height = 500)

stat.test <- compare_means(
  c(X1, X2) ~ cluster, data = df,
  method = "t.test"
)
stat.test <- stat.test |>
  mutate(y.position = c(5, 5))

# Add manually p-values from stat.test data
# First specify the y.position of each comparison
# p + stat_pvalue_manual(stat.test, label = "p.adj")

# Create box plot using ggplot2
p1 <- ggboxplot(df, x = "cluster", y = "X1",
                color = "cluster", palette = get_palette("Accent", 3)) +
  stat_pvalue_manual(stat.test, label = "p.format")
p2 <- ggboxplot(df, x = "cluster", y = "X2",
                color = "cluster", palette = get_palette("Accent", 3)) +
  stat_pvalue_manual(stat.test, label = "p.format")

ggarrange(p1, p2, common.legend = TRUE) |>
  ggexport(filename = "../seminar_paper-custom/figures/fig-boxplot-S15.png",
           width = 800, height = 400)

# ggarrange(p_mvn_gray, p_mvn_clustered, p1, p2)
