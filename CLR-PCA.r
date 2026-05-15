# Install packages
my_packages <- c()

install.packages(my_packages)

# Check for missing files and install them
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(RColorBrewer)
library(patchwork)

#Load Files
metaphlan <- read.table(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header       = TRUE,
  sep          = "\t",
  row.names    = 1,
  check.names  = FALSE,
)

metadata <- read.table(
  "metadata_merged.txt", #This is a placeholder for your merged metaphlan file
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
colnames(metadata)[1] <- "SampleID"
 
cat("Metadata loaded:", nrow(metadata), "samples x",
    ncol(metadata), "columns\n")

# Species-level only
metaphlan <- metaphlan[grepl("s__", rownames(metaphlan)), ]
 
# Transpose: samples as rows, taxa as columns
abund_mat <- as.data.frame(t(metaphlan))
 
cat("MetaPhlAn table loaded:", nrow(abund_mat), "samples x",
    ncol(abund_mat), "taxa\n")

# Align Samples
common_samples  <- intersect(rownames(abund_mat), metadata$SampleID)
cat("Common samples:", length(common_samples), "\n")
 
abund_mat       <- abund_mat[common_samples, ]
metadata <- metadata[match(common_samples,
                                         metadata$SampleID), ]
rownames(metadata) <- metadata$SampleID

# Prevalence Filter
prevalence_threshold <- 0.20
prev      <- colMeans(abund_mat > 0)
taxa_keep <- names(prev[prev >= prevalence_threshold])
cat("Taxa kept after prevalence filter:", length(taxa_keep), "/",
    ncol(abund_mat), "\n")
 
abund_filt <- abund_mat[, taxa_keep]

# CLR Transformation
clr_transform <- function(x) {
  x <- x + 1e-6
  log(x) - mean(log(x))
}
 
abund_clr <- as.data.frame(t(apply(abund_filt, 1, clr_transform)))

# PCA
pca_result <- prcomp(abund_clr, center = TRUE, scale. = FALSE)
 
# Variance explained
var_explained <- summary(pca_result)$importance
pct_var       <- round(var_explained[2, ] * 100, 1)
 
cat("\nVariance explained:\n")
cat("  PC1:", pct_var[1], "%\n")
cat("  PC2:", pct_var[2], "%\n")
cat("  PC3:", pct_var[3], "%\n")
 
# PCA scores
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleID <- rownames(pca_scores)
 
# Merge with metadata
pca_df <- left_join(pca_scores, metadata, by = "SampleID")

# Set Grouping Variable
group_col <- "Group"
 
if (!group_col %in% colnames(pca_df)) {
  cat("Warning: '", group_col, "' not found. Using no grouping.\n", sep = "")
  pca_df$Group <- "All"
  group_col    <- "Group"
}
 
pca_df$GroupVar <- pca_df[[group_col]]

# Plots
n_groups <- length(unique(pca_df$GroupVar))
pal      <- if (n_groups <= 8) brewer.pal(max(3, n_groups), "Set1") else
              colorRampPalette(brewer.pal(8, "Set1"))(n_groups)

# PC1 v.s. PC2
p1 <- ggplot(pca_df, aes(x = PC1, y = PC2,
                          color = GroupVar, fill = GroupVar)) +
  stat_ellipse(aes(group = GroupVar),
               geom = "polygon", alpha = 0.10,
               linetype = "dashed", level = 0.95) +
  geom_point(size = 3.5, alpha = 0.85, shape = 21,
             color = "white", stroke = 0.5,
             aes(fill = GroupVar)) +
  geom_text_repel(aes(label = SampleID),
                  size = 2.5, color = "grey30",
                  max.overlaps = 15, show.legend = FALSE) +
  scale_color_manual(values = pal, name = group_col) +
  scale_fill_manual(values  = pal, name = group_col) +
  labs(
    title    = "CLR-PCA — PC1 vs PC2",
    subtitle = paste0("PC1: ", pct_var[1], "%  |  PC2: ", pct_var[2], "%"),
    x        = paste0("PC1 (", pct_var[1], "%)"),
    y        = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# PC1 v.s. PC3 
p2 <- ggplot(pca_df, aes(x = PC1, y = PC3,
                          color = GroupVar, fill = GroupVar)) +
  stat_ellipse(aes(group = GroupVar),
               geom = "polygon", alpha = 0.10,
               linetype = "dashed", level = 0.95) +
  geom_point(size = 3.5, alpha = 0.85, shape = 21,
             color = "white", stroke = 0.5,
             aes(fill = GroupVar)) +
  scale_color_manual(values = pal, name = group_col) +
  scale_fill_manual(values  = pal, name = group_col) +
  labs(
    title = "CLR-PCA — PC1 vs PC3",
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC3 (", pct_var[3], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# Scree Plot
scree_df <- data.frame(
  PC  = paste0("PC", 1:min(15, length(pct_var))),
  Var = pct_var[1:min(15, length(pct_var))]
)
scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)
 
p3 <- ggplot(scree_df, aes(x = PC, y = Var)) +
  geom_col(fill = "#4393C3", alpha = 0.85, width = 0.7) +
  geom_line(aes(group = 1), color = "#D6604D", linewidth = 0.8) +
  geom_point(color = "#D6604D", size = 2.5) +
  labs(
    title = "Scree Plot",
    x     = "Principal Component",
    y     = "Variance Explained (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

# Top Taxa Loadings (PC1 & PC2)
loadings <- as.data.frame(pca_result$rotation[, 1:4])
loadings$Taxon <- rownames(loadings)
 
# Shorten taxon names for readability (keep species epithet)
loadings$TaxonShort <- gsub(".*s__", "", loadings$Taxon)
 
top_n <- 15
 
top_pc1 <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  slice_head(n = top_n) %>%
  mutate(TaxonShort = fct_reorder(TaxonShort, PC1))
 
top_pc2 <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n = top_n) %>%
  mutate(TaxonShort = fct_reorder(TaxonShort, PC2))
 
p4 <- ggplot(top_pc1, aes(x = PC1, y = TaxonShort,
                            fill = PC1 > 0)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#D6604D",
                                "FALSE" = "#2166AC"),
                    guide = "none") +
  labs(
    title = paste0("Top ", top_n, " Taxa — PC1 Loadings"),
    x     = paste0("PC1 loading (", pct_var[1], "%)"),
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )
 
p5 <- ggplot(top_pc2, aes(x = PC2, y = TaxonShort,
                            fill = PC2 > 0)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#D6604D",
                                "FALSE" = "#2166AC"),
                    guide = "none") +
  labs(
    title = paste0("Top ", top_n, " Taxa — PC2 Loadings"),
    x     = paste0("PC2 loading (", pct_var[2], "%)"),
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

# Biplot (Samples & Top Taxon Arrows)
top_biplot <- loadings %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(magnitude)) %>%
  slice_head(n = 10)
 
# Scale arrows to sample score range
scale_factor <- max(abs(pca_df[, c("PC1", "PC2")])) /
                max(top_biplot$magnitude) * 0.6
 
p6 <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = GroupVar), size = 3, shape = 21,
             color = "white", stroke = 0.5, alpha = 0.8) +
  geom_segment(
    data = top_biplot,
    aes(x = 0, y = 0,
        xend = PC1 * scale_factor,
        yend = PC2 * scale_factor),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "grey30", linewidth = 0.5
  ) +
  geom_text_repel(
    data = top_biplot,
    aes(x = PC1 * scale_factor,
        y = PC2 * scale_factor,
        label = TaxonShort),
    size = 2.5, color = "grey20", fontface = "italic",
    max.overlaps = 20
  ) +
  scale_fill_manual(values = pal, name = group_col) +
  labs(
    title    = "CLR-PCA Biplot",
    subtitle = "Top 10 taxa by loading magnitude",
    x        = paste0("PC1 (", pct_var[1], "%)"),
    y        = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey50"),
    panel.grid.minor = element_blank()
  )

# Save all plots
# Combined overview (PC1v2, PC1v3, scree)
pdf("clr_pca_overview.pdf", width = 16, height = 6)
print(p1 + p2 + p3)
dev.off()
cat("Saved: clr_pca_overview.pdf\n")
 
# Loadings
pdf("clr_pca_loadings.pdf", width = 14, height = 6)
print(p4 + p5)
dev.off()
cat("Saved: clr_pca_loadings.pdf\n")
 
# Biplot
ggsave("clr_pca_biplot.pdf", p6, width = 10, height = 8, dpi = 300)
cat("Saved: clr_pca_biplot.pdf\n")

# PERMANOVA (group effect on community composition)
if (group_col %in% colnames(metadata_merged)) {
  set.seed(42)
  perm_result <- adonis2(
    abund_clr ~ metadata_merged[[group_col]],
    method  = "euclidean",
    permutations = 999
  )
  print(perm_result)
 
  # Save PERMANOVA result
  write.csv(as.data.frame(perm_result),
            "clr_pca_permanova.csv", row.names = TRUE)
  cat("Saved: clr_pca_permanova.csv\n")
}

#Save PCA Scores and Loadings
write.csv(pca_df,    "clr_pca_scores.csv",   row.names = FALSE)
write.csv(loadings,  "clr_pca_loadings.csv", row.names = FALSE)
cat("\nSaved: clr_pca_scores.csv, clr_pca_loadings.csv\n")
