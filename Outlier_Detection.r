#This code detects both sample and taxon level outliers using PCA-based (Mahalanobis distance) method. 
# Install Packages
my_packages <- c("tidyverse", "vegan", "ggplot2", "ggrepel",
                   "RColorBrewer", "patchwork", "factoextra", "mvoutlier")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(factoextra)
library(mvoutlier)

#Load Files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE, 
)

metadata <- read.delim(
  "metadata_merged.txt",  #This is a placeholder for your merged metadata file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE
)
 
meta_id_col <- colnames(metadata)[1]

# Filter to species level & build abundance matrix
species_data <- metaphlan %>%
  filter(grepl("s__", .[[1]]) & !grepl("t__", .[[1]]))
 
colnames(species_data)[1] <- "taxonomy"
species_data <- species_data %>%
  mutate(species = gsub(".*\\|", "", taxonomy)) %>%
  select(species, everything(), -taxonomy)

# Matrix: rows = samples, cols = species
abund_mat <- species_data %>%
  column_to_rownames("species") %>%
  t() %>%
  as.data.frame()
 
# Remove species present in <10 % of samples (reduce noise)
prev_filter <- colSums(abund_mat > 0) / nrow(abund_mat) >= 0.10
abund_mat   <- abund_mat[, prev_filter]
message(sprintf("Retained %d species (prevalence ≥ 10%%) across %d samples.",
                ncol(abund_mat), nrow(abund_mat)))
 
# Hellinger-transform before PCA (standard for compositional data)
abund_hell <- decostand(abund_mat, method = "hellinger")

# PCA
pca_res  <- prcomp(abund_hell, center = TRUE, scale. = FALSE)
pca_df   <- as.data.frame(pca_res$x)
pca_df$sample_id <- rownames(pca_df)
 
# Variance explained
var_exp  <- summary(pca_res)$importance[2, 1:2] * 100   # PC1, PC2
 
# ── 5. SAMPLE-LEVEL OUTLIERS via Mahalanobis Distance ────────
# Use top PCs that explain ≥ 80 % cumulative variance
cum_var   <- cumsum(summary(pca_res)$importance[2, ])
n_pc      <- max(2, which(cum_var >= 0.80)[1])
pc_scores <- pca_res$x[, 1:n_pc]

# cov_mat       <- cov(pc_scores)
center_vec    <- colMeans(pc_scores)
maha_dist     <- mahalanobis(pc_scores, center = center_vec, cov = cov_mat)
maha_p        <- pchisq(maha_dist, df = n_pc, lower.tail = FALSE)
maha_p_adj    <- p.adjust(maha_p, method = "BH")
 
sample_outliers <- data.frame(
  sample_id      = rownames(pc_scores),
  mahalanobis_D2 = maha_dist,
  p_value        = maha_p,
  p_adj_BH       = maha_p_adj,
  outlier        = maha_p_adj < 0.05
)
 
n_out <- sum(sample_outliers$outlier)
message(sprintf("Sample outliers detected: %d", n_out))

# PCA Coloured by Mahalanobis Distance Plot
plot_df <- pca_df %>%
  left_join(sample_outliers, by = "sample_id") %>%
  left_join(metadata, by = setNames(meta_id_col, "sample_id"))
 
# Auto-detect grouping variable
group_col <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) & n_distinct(.) > 1 & n_distinct(.) < 20)) %>%
  colnames() %>% first()
 
base_pca <- ggplot(plot_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = mahalanobis_D2, shape = outlier),
             size = 4, stroke = 0.6, color = "grey30") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
    name    = "Mahalanobis D²"
  ) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24),
                     labels = c("Normal", "Outlier"),
                     name   = NULL) +
  geom_text_repel(
    data   = filter(plot_df, outlier),
    aes(label = sample_id), size = 3, color = "red3",
    box.padding = 0.5, max.overlaps = 20
  ) +
  labs(
    title    = "PCA – Sample-Level Outlier Detection",
    subtitle = sprintf("Mahalanobis D²  |  %d PCs  |  Outliers (BH-adj p < 0.05): %d", n_pc, n_out),
    x        = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp[2])
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
 
# Overlay group ellipses if group variable found
if (!is.null(group_col) && !is.na(group_col)) {
  n_grp <- n_distinct(plot_df[[group_col]])
  pal   <- colorRampPalette(brewer.pal(min(n_grp, 8), "Dark2"))(n_grp)
  base_pca <- base_pca +
    stat_ellipse(aes(color = .data[[group_col]]), linewidth = 0.8,
                 level = 0.95, linetype = "dashed") +
    scale_color_manual(values = pal, name = group_col)
}
 
# Mahalanobis distance distribution plot
p_maha_hist <- ggplot(sample_outliers, aes(x = mahalanobis_D2, fill = outlier)) +
  geom_histogram(bins = 30, color = "white", alpha = 0.85) +
  scale_fill_manual(values = c(`FALSE` = "#4E9AF1", `TRUE` = "#E53935"),
                    labels = c("Normal", "Outlier")) +
  geom_vline(
    xintercept = qchisq(0.95, df = n_pc),
    linetype = "dashed", color = "grey30", linewidth = 0.8
  ) +
  annotate("text", x = qchisq(0.95, df = n_pc), y = Inf,
           label = " χ² 95th pctl", hjust = 0, vjust = 2, size = 3.2) +
  labs(title = "Mahalanobis D² Distribution",
       x = "Mahalanobis D²", y = "Count", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
 
sample_panel <- base_pca / p_maha_hist + plot_layout(heights = c(2, 1))
ggsave("outlier_sample_level.pdf", sample_panel, width = 11, height = 13, dpi = 300)
message("Saved: outlier_sample_level.pdf")

# AXON-LEVEL OUTLIERS via Mahalanobis Distance
# Transpose: rows = species, cols = samples (PC space)
abund_t     <- t(abund_hell)                      # species × samples
species_pca <- prcomp(abund_t, center = TRUE, scale. = FALSE)
 
cum_var_sp  <- cumsum(summary(species_pca)$importance[2, ])
n_pc_sp     <- max(2, which(cum_var_sp >= 0.80)[1])
sp_scores   <- species_pca$x[, 1:n_pc_sp]
 
var_exp_sp  <- summary(species_pca)$importance[2, 1:2] * 100
 
cov_sp      <- cov(sp_scores)
center_sp   <- colMeans(sp_scores)
maha_sp     <- mahalanobis(sp_scores, center = center_sp, cov = cov_sp)
p_sp        <- pchisq(maha_sp, df = n_pc_sp, lower.tail = FALSE)
p_sp_adj    <- p.adjust(p_sp, method = "BH")
 
taxon_outliers <- data.frame(
  species        = rownames(sp_scores),
  mahalanobis_D2 = maha_sp,
  p_value        = p_sp,
  p_adj_BH       = p_sp_adj,
  mean_abundance = rowMeans(abund_t),
  outlier        = p_sp_adj < 0.05
)
 
n_tax_out <- sum(taxon_outliers$outlier)
message(sprintf("Taxon outliers detected: %d", n_tax_out))
 
# Taxon PCA plot
taxon_pca_df <- as.data.frame(sp_scores) %>%
  rownames_to_column("species") %>%
  left_join(taxon_outliers %>% select(species, mahalanobis_D2, outlier), by = "species")
 
p_taxon_pca <- ggplot(taxon_pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = mahalanobis_D2, shape = outlier),
             size = 3.5, stroke = 0.5, color = "grey30") +
  scale_fill_gradientn(
    colours = c("#1B7837", "#7FBF7B", "#F7F7F7", "#AF8DC3", "#762A83"),
    name    = "Mahalanobis D²"
  ) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24),
                     labels = c("Normal", "Outlier"), name = NULL) +
  geom_text_repel(
    data   = filter(taxon_pca_df, outlier),
    aes(label = gsub("s__", "", species)), size = 2.8, color = "purple4",
    box.padding = 0.5, max.overlaps = 25
  ) +
  labs(
    title    = "PCA – Taxon-Level Outlier Detection",
    subtitle = sprintf("Mahalanobis D²  |  %d PCs  |  Outlier taxa (BH-adj p < 0.05): %d", n_pc_sp, n_tax_out),
    x        = sprintf("PC1 (%.1f%%)", var_exp_sp[1]),
    y        = sprintf("PC2 (%.1f%%)", var_exp_sp[2])
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
 
# Top outlier taxa bar plot
p_taxon_bar <- taxon_outliers %>%
  slice_max(mahalanobis_D2, n = 20) %>%
  mutate(species = gsub("s__", "", species),
         species = fct_reorder(species, mahalanobis_D2)) %>%
  ggplot(aes(x = mahalanobis_D2, y = species, fill = outlier)) +
  geom_col(color = "white", alpha = 0.85) +
  scale_fill_manual(values = c(`FALSE` = "#7986CB", `TRUE` = "#E53935"),
                    labels = c("Normal", "Outlier"), name = NULL) +
  geom_vline(xintercept = qchisq(0.95, df = n_pc_sp),
             linetype = "dashed", color = "grey30") +
  labs(title = "Top 20 Taxa by Mahalanobis D²",
       x = "Mahalanobis D²", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
 
taxon_panel <- p_taxon_pca / p_taxon_bar + plot_layout(heights = c(2, 1.5))
ggsave("outlier_taxon_level.pdf", taxon_panel, width = 11, height = 14, dpi = 300)
message("Saved: outlier_taxon_level.pdf")

# Combined Summary Plot
combined <- (base_pca | p_taxon_pca) /
            (p_maha_hist | p_taxon_bar) +
  plot_annotation(
    title    = "Microbiome Outlier Detection – PCA + Mahalanobis Distance",
    subtitle = sprintf("Samples: %d  |  Taxa: %d  |  Sample outliers: %d  |  Taxon outliers: %d",
                       nrow(abund_mat), ncol(abund_mat), n_out, n_tax_out),
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10))
  )
 
ggsave("outlier_combined.pdf", combined, width = 16, height = 14, dpi = 300)
message("Saved: outlier_combined.pdf")

# Save Tables
write.csv(sample_outliers %>% arrange(p_adj_BH),
          "outlier_samples.csv",  row.names = FALSE)
write.csv(taxon_outliers  %>% arrange(p_adj_BH),
          "outlier_taxa.csv",     row.names = FALSE)

message("\n✓ Outlier detection complete.")
