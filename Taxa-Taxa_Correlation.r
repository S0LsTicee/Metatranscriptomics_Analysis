#Install packages
my_packages <- c("tidyverse", "Hmisc", "corrplot", "ggplot2",
              "ggrepel", "RColorBrewer", "patchwork", "igraph",
              "ggraph", "reshape2", "scales", "pheatmap")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(Hmisc)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(igraph)
library(ggraph)
library(reshape2)
library(scales)
library(pheatmap)

#Load Files
metaphlan <- read.table(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header       = TRUE,
  sep          = "\t",
  row.names    = 1,
  check.names  = FALSE,
)

metadata <- read.table(
  "metadata_merged.txt", #This is a placeholder for your merged metadata file
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
colnames(metadata)[1] <- "SampleID"
rownames(metadata)    <- metadata$SampleID

# Species-level only
metaphlan <- metaphlan[grepl("s__", rownames(metaphlan)), ]
rownames(metaphlan) <- gsub(".*s__", "", rownames(metaphlan))

# Align Samples
common_samples  <- intersect(colnames(metaphlan),
                              metadata$SampleID)
cat("Common samples:", length(common_samples), "\n")
 
abund_mat       <- metaphlan[, common_samples]   # taxa x samples
metadata <- metadata[common_samples, ]

# Set Grouping Variable
group_col <- "Group"
 
if (!group_col %in% colnames(metadata)) {
  cat("Warning: '", group_col, "' not found — using 'All'.\n", sep = "")
  metadata$Group <- "All"
  group_col <- "Group"
}
 
groups     <- unique(metadata[[group_col]])
n_groups   <- length(groups)
pal        <- if (n_groups <= 8) brewer.pal(max(3, n_groups), "Set1") else
                colorRampPalette(brewer.pal(8, "Set1"))(n_groups)
names(pal) <- groups

# Prevalence filter
prev_threshold <- 0.20
prev           <- rowMeans(abund_mat > 0)
taxa_keep      <- names(prev[prev >= prev_threshold])
cat("Taxa after prevalence filter:", length(taxa_keep), "/",
    nrow(abund_mat), "\n")
 
abund_filt <- abund_mat[taxa_keep, ]

# CLR transformation (samples x taxa)
clr_transform <- function(x) {
  x <- x + 1e-6
  log(x) - mean(log(x))
}
abund_clr <- as.data.frame(t(apply(t(abund_filt), 1, clr_transform)))

# Raw relative abundance (samples x taxa), for Spearman/Pearson
abund_ra <- as.data.frame(t(abund_filt))   # samples x taxa

# Correlation Methods
# Helper: compute correlation + BH-corrected p-values
compute_cor <- function(mat, method = "spearman") {
  res      <- rcorr(as.matrix(mat), type = method)
  r        <- res$r
  p        <- res$P
  pv       <- p[upper.tri(p)]
  pa       <- p.adjust(pv, method = "BH")
  padj     <- matrix(1, nrow(p), ncol(p))
  padj[upper.tri(padj)] <- pa
  padj     <- padj + t(padj)
  diag(padj) <- 1
  rownames(padj) <- colnames(padj) <- colnames(mat)
  list(r = r, p = p, padj = padj)
}

# Spearman on Raw Relative Abundance
spearman_res <- compute_cor(abund_ra, method = "spearman")

# Pearson on CLR-transformed Data
pearson_clr_res <- compute_cor(abund_clr, method = "pearson")

# Combined Results Table
taxa <- colnames(abund_clr)
n    <- length(taxa)
 
pairs_df <- do.call(rbind, lapply(1:(n - 1), function(i) {
  lapply((i + 1):n, function(j) {
    data.frame(
      Taxon1          = taxa[i],
      Taxon2          = taxa[j],
      Spearman_r      = spearman_res$r[i, j],
      Spearman_padj   = spearman_res$padj[i, j],
      Pearson_CLR_r   = pearson_clr_res$r[i, j],
      Pearson_CLR_padj= pearson_clr_res$padj[i, j],
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
})) %>% bind_rows()
 
# Classify significance
pairs_df <- pairs_df %>%
  mutate(
    Sig_Spearman   = Spearman_padj   < 0.05,
    Sig_Pearson_CLR= Pearson_CLR_padj < 0.05,
    Sig_Both       = Sig_Spearman & Sig_Pearson_CLR,
    Direction      = case_when(
      Sig_Both & Spearman_r > 0 ~ "Positive (both)",
      Sig_Both & Spearman_r < 0 ~ "Negative (both)",
      Sig_Spearman              ~ "Spearman only",
      Sig_Pearson_CLR           ~ "CLR-Pearson only",
      TRUE                      ~ "Not significant"
    )
  ) %>%
  arrange(Spearman_padj)

cat("Top 10 correlated pairs (Spearman):\n")
print(head(pairs_df %>%
             filter(Sig_Spearman) %>%
             select(Taxon1, Taxon2, Spearman_r, Spearman_padj,
                    Pearson_CLR_r, Pearson_CLR_padj), 10))

# Group-stratified Correlations
group_cor_list <- lapply(groups, function(grp) {
  samps <- rownames(metadata[metadata[[group_col]] == grp, ])
  sub   <- abund_clr[samps, ]
  res   <- compute_cor(sub, method = "pearson")
  res$r
})
names(group_cor_list) <- groups
 
# Two-group Differential Correlation
if (n_groups == 2) {
  diff_cor_mat <- group_cor_list[[2]] - group_cor_list[[1]]
  cat("Differential correlation matrix computed (",
      groups[2], "-", groups[1], ").\n")
}

# Spearman Heatmap
pdf("taxa_cor_spearman_heatmap.pdf", width = 13, height = 12)
corrplot(
  spearman_res$r,
  method    = "color",
  type      = "upper",
  order     = "hclust",
  addrect   = 4,
  col       = colorRampPalette(c("#2166AC", "white", "#D6604D"))(200),
  tl.cex    = 0.6,
  tl.col    = "black",
  p.mat     = spearman_res$padj,
  sig.level = 0.05,
  insig     = "blank",
  title     = "Taxa-Taxa Spearman Correlation (FDR < 0.05)",
  mar       = c(0, 0, 2, 0)
)
dev.off()
cat("Saved: taxa_cor_spearman_heatmap.pdf\n")

# CLR-Pearson Heatmap
pdf("taxa_cor_pearson_clr_heatmap.pdf", width = 13, height = 12)
corrplot(
  pearson_clr_res$r,
  method    = "color",
  type      = "upper",
  order     = "hclust",
  addrect   = 4,
  col       = colorRampPalette(c("#2166AC", "white", "#D6604D"))(200),
  tl.cex    = 0.6,
  tl.col    = "black",
  p.mat     = pearson_clr_res$padj,
  sig.level = 0.05,
  insig     = "blank",
  title     = "Taxa-Taxa Pearson Correlation on CLR (FDR < 0.05)",
  mar       = c(0, 0, 2, 0)
)
dev.off()
cat("Saved: taxa_cor_pearson_clr_heatmap.pdf\n")

# Method Comparison Scatter
p_compare <- ggplot(pairs_df,
                    aes(x = Spearman_r, y = Pearson_CLR_r,
                        color = Direction)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  geom_smooth(data = pairs_df, aes(x = Spearman_r, y = Pearson_CLR_r),
              method = "lm", se = TRUE, color = "grey30",
              inherit.aes = FALSE, linewidth = 0.7) +
  scale_color_manual(
    values = c(
      "Positive (both)"    = "#D73027",
      "Negative (both)"    = "#4575B4",
      "Spearman only"      = "#FC8D59",
      "CLR-Pearson only"   = "#74ADD1",
      "Not significant"    = "grey80"
    ),
    name = "Significance"
  ) +
  labs(
    title    = "Spearman vs CLR-Pearson Correlation",
    subtitle = paste0("All ", nrow(pairs_df), " taxa pairs"),
    x        = "Spearman r (raw RA)",
    y        = "Pearson r (CLR-transformed)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Correlation Distributions
dist_df <- bind_rows(
  data.frame(r = pairs_df$Spearman_r,    Method = "Spearman (raw RA)"),
  data.frame(r = pairs_df$Pearson_CLR_r, Method = "Pearson (CLR)")
)
 
p_dist <- ggplot(dist_df, aes(x = r, fill = Method, color = Method)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values  = c("Spearman (raw RA)" = "#D6604D",
                                 "Pearson (CLR)"     = "#4393C3")) +
  scale_color_manual(values = c("Spearman (raw RA)" = "#D6604D",
                                 "Pearson (CLR)"     = "#4393C3")) +
  labs(
    title = "Distribution of Pairwise Correlation Coefficients",
    x     = "Correlation coefficient (r)",
    y     = "Density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.title     = element_blank()
  )

# Spearman Top Positive & Negative Pairs
top_pos <- pairs_df %>% filter(Sig_Spearman, Spearman_r > 0) %>%
  arrange(desc(Spearman_r)) %>% slice_head(n = 15) %>%
  mutate(Pair = paste(Taxon1, "—", Taxon2),
         Pair = fct_reorder(Pair, Spearman_r))
 
top_neg <- pairs_df %>% filter(Sig_Spearman, Spearman_r < 0) %>%
  arrange(Spearman_r) %>% slice_head(n = 15) %>%
  mutate(Pair = paste(Taxon1, "—", Taxon2),
         Pair = fct_reorder(Pair, Spearman_r))
 
p_top_pos <- ggplot(top_pos, aes(x = Spearman_r, y = Pair,
                                   fill = Spearman_padj)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_gradient(low = "#D73027", high = "#FC8D59",
                      name = "FDR", trans = "log10") +
  labs(title = "Top Positive Pairs (Spearman)",
       x = "Spearman r", y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 7.5),
        panel.grid.major.y = element_blank())
 
p_top_neg <- ggplot(top_neg, aes(x = Spearman_r, y = Pair,
                                   fill = Spearman_padj)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_gradient(low = "#4575B4", high = "#74ADD1",
                      name = "FDR", trans = "log10") +
  labs(title = "Top Negative Pairs (Spearman)",
       x = "Spearman r", y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 7.5),
        panel.grid.major.y = element_blank())

# Pheatmap of Significant Taxa 
sig_taxa <- unique(c(
  pairs_df$Taxon1[pairs_df$Sig_Both],
  pairs_df$Taxon2[pairs_df$Sig_Both]
))
 
if (length(sig_taxa) >= 4) {
  sub_mat <- pearson_clr_res$r[sig_taxa, sig_taxa]
 
  # Annotation: mean abundance per taxon
  ann_row <- data.frame(
    MeanAbund = rowMeans(abund_filt[sig_taxa, ]),
    row.names = sig_taxa
  )
 
  pdf("taxa_cor_sig_pheatmap.pdf", width = 12, height = 11)
  pheatmap(
    sub_mat,
    color            = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
    breaks           = seq(-1, 1, length.out = 101),
    clustering_method = "ward.D2",
    annotation_row   = ann_row,
    fontsize_row     = 7,
    fontsize_col     = 7,
    border_color     = NA,
    main             = "CLR-Pearson Correlations — Significant Taxa Only"
  )
  dev.off()
  cat("Saved: taxa_cor_sig_pheatmap.pdf\n")
}

# Correlation Network (CLR-Pearson, FDR < 0.05)
net_pairs <- pairs_df %>%
  filter(Sig_Both, abs(Pearson_CLR_r) >= 0.40)
 
if (nrow(net_pairs) > 0) {
  g <- graph_from_data_frame(
    data.frame(from   = net_pairs$Taxon1,
               to     = net_pairs$Taxon2,
               weight = net_pairs$Pearson_CLR_r),
    directed = FALSE,
    vertices = taxa
  )
  E(g)$sign       <- ifelse(E(g)$weight > 0, "positive", "negative")
  E(g)$abs_weight <- abs(E(g)$weight)
  V(g)$degree     <- degree(g)
  V(g)$betweenness <- betweenness(g, normalized = TRUE)
  g_sub <- delete_vertices(g, V(g)[degree(g) == 0])
 
  set.seed(42)
  p_net <- ggraph(g_sub, layout = "fr") +
    geom_edge_link(aes(color = sign, width = abs_weight,
                       alpha = abs_weight)) +
    scale_edge_color_manual(
      values = c("positive" = "#D6604D", "negative" = "#2166AC"),
      name   = "Direction"
    ) +
    scale_edge_width(range = c(0.3, 2.5), guide = "none") +
    scale_edge_alpha(range = c(0.4, 0.9), guide = "none") +
    geom_node_point(aes(size = degree, fill = betweenness),
                    shape = 21, color = "white", stroke = 0.5) +
    scale_size(range = c(2.5, 9), name = "Degree") +
    scale_fill_gradient(low = "#FEE08B", high = "#1A9850",
                        name = "Betweenness\n(normalized)") +
    geom_node_text(aes(label = name), repel = TRUE,
                   size = 2.4, color = "grey20", fontface = "italic") +
    theme_graph() +
    labs(
      title    = "Taxa-Taxa Correlation Network",
      subtitle = paste0("CLR-Pearson |r| ≥ 0.40  &  FDR < 0.05 (both methods)  |  ",
                        ecount(g_sub), " edges")
    ) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey50"))
 
  ggsave("taxa_cor_network.pdf", p_net, width = 14, height = 10, dpi = 300)
  cat("Saved: taxa_cor_network.pdf\n")
}

# Two-Groups Differential Correlation Heatmaps
if (n_groups == 2 && exists("diff_cor_mat")) {
  pdf("taxa_cor_differential.pdf", width = 13, height = 12)
  corrplot(
    diff_cor_mat,
    method  = "color",
    type    = "upper",
    order   = "hclust",
    col     = colorRampPalette(c("#762A83", "white", "#1B7837"))(200),
    tl.cex  = 0.6,
    tl.col  = "black",
    title   = paste0("Differential Correlation (CLR-Pearson): ",
                      groups[2], " − ", groups[1]),
    mar     = c(0, 0, 2, 0),
    cl.lim  = c(-1, 1)
  )
  dev.off()
  cat("Saved: taxa_cor_differential.pdf\n")
}

# Save Plots & Tables
pdf("taxa_cor_method_comparison.pdf", width = 14, height = 6)
print(p_compare + p_dist)
dev.off()
cat("Saved: taxa_cor_method_comparison.pdf\n")
 
pdf("taxa_cor_top_pairs.pdf", width = 14, height = 7)
print(p_top_pos + p_top_neg)
dev.off()
cat("Saved: taxa_cor_top_pairs.pdf\n")
 
write.csv(pairs_df, "taxa_cor_all_pairs.csv",         row.names = FALSE)
write.csv(
  pairs_df %>% filter(Sig_Both),
  "taxa_cor_sig_pairs.csv", row.names = FALSE
)
write.csv(
  as.data.frame(pearson_clr_res$r),
  "taxa_cor_clr_pearson_matrix.csv", row.names = TRUE
)
write.csv(
  as.data.frame(spearman_res$r),
  "taxa_cor_spearman_matrix.csv", row.names = TRUE
)
 
cat("\nSaved tables:\n")
