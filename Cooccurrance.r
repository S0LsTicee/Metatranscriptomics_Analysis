#Install packages
my_packages <- c("tidyverse", "vegan", "igraph", "ggraph",
              "corrplot", "Hmisc", "reshape2", "RColorBrewer")

install.packages(my_packages)

# Check for missing files and install them
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(vegan)
library(igraph)
library(ggraph)
library(corrplot)
library(Hmisc)
library(reshape2)
library(RColorBrewer)

# Load Files
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
 
cat("Metadata loaded:", nrow(metadata), "samples x",
    ncol(metadata), "columns\n")

# Keep only species-level rows (containing "s__")
metaphlan <- metaphlan[grepl("s__", rownames(metaphlan)), ]
 
# Transpose: samples as rows, taxa as columns
abund_mat <- as.data.frame(t(metaphlan))
 
cat("MetaPhlAn table loaded:", nrow(abund_mat), "samples x",
    ncol(abund_mat), "taxa\n")

# Align Samples, assuming the first column of metadata is the sample ID
colnames(metadata)[1] <- "SampleID"
 
common_samples  <- intersect(rownames(abund_mat), metadata$SampleID)
cat("Common samples:", length(common_samples), "\n")
 
abund_mat       <- abund_mat[common_samples, ]
metadata <- metadata[match(common_samples,
                                         metadata$SampleID), ]

# Prevalence filter: keep taxa present in >= 20% of samples
prevalence_threshold <- 0.20
prev      <- colMeans(abund_mat > 0)
taxa_keep <- names(prev[prev >= prevalence_threshold])
cat("Taxa kept after prevalence filter:", length(taxa_keep), "/",
    ncol(abund_mat), "\n")
 
abund_filt <- abund_mat[, taxa_keep]

# CLR transformation (compositional-data aware)
clr_transform <- function(x) {
  x <- x + 1e-6          # pseudocount for zeros
  log(x) - mean(log(x))
}
 
abund_clr <- as.data.frame(t(apply(abund_filt, 1, clr_transform)))

# Spearman Correlation & FDR Correction
cor_result <- rcorr(as.matrix(abund_clr), type = "spearman")
cor_mat    <- cor_result$r
pval_mat   <- cor_result$P
 
## BH FDR correction
pval_vec <- pval_mat[upper.tri(pval_mat)]
padj_vec <- p.adjust(pval_vec, method = "BH")
padj_mat <- matrix(0, nrow = nrow(pval_mat), ncol = ncol(pval_mat))
padj_mat[upper.tri(padj_mat)] <- padj_vec
padj_mat <- padj_mat + t(padj_mat)
diag(padj_mat) <- 1
rownames(padj_mat) <- colnames(padj_mat) <- colnames(abund_clr)
 
cat("Correlation matrix computed for", ncol(abund_clr), "taxa.\n")

# Correlation Heatmap
pdf("cooccurrence_heatmap.pdf", width = 12, height = 11)
 
corrplot(
  cor_mat,
  method    = "color",
  type      = "upper",
  order     = "hclust",
  addrect   = 4,
  col       = colorRampPalette(c("#2166AC", "white", "#D6604D"))(200),
  tl.cex    = 0.6,
  tl.col    = "black",
  p.mat     = padj_mat,
  sig.level = 0.05,
  insig     = "blank",
  title     = "Microbial Co-occurrence (Spearman, FDR < 0.05)",
  mar       = c(0, 0, 2, 0)
)
 
dev.off()
cat("Saved: cooccurrence_heatmap.pdf\n")

# Build Cooccurrance Network
r_threshold   <- 0.40
fdr_threshold <- 0.05
 
edge_idx <- which(
  abs(cor_mat) >= r_threshold &
    padj_mat   <  fdr_threshold &
    upper.tri(cor_mat),
  arr.ind = TRUE
)
 
edges <- data.frame(
  from   = rownames(cor_mat)[edge_idx[, 1]],
  to     = colnames(cor_mat)[edge_idx[, 2]],
  weight = cor_mat[edge_idx],
  stringsAsFactors = FALSE
)
 
cat("Edges in network (|r| >=", r_threshold,
    "& FDR <", fdr_threshold, "):", nrow(edges), "\n")
 
g <- graph_from_data_frame(edges, directed = FALSE,
                            vertices = colnames(abund_clr))
 
E(g)$sign        <- ifelse(E(g)$weight > 0, "positive", "negative")
E(g)$abs_weight  <- abs(E(g)$weight)
V(g)$degree      <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)
 
node_stats <- data.frame(
  Taxon       = V(g)$name,
  Degree      = V(g)$degree,
  Betweenness = round(V(g)$betweenness, 4)
) %>% arrange(desc(Degree))

# Network Plot
set.seed(123)
 
p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(color = sign, width = abs_weight, alpha = abs_weight)
  ) +
  scale_edge_color_manual(
    values = c("positive" = "#D6604D", "negative" = "#2166AC"),
    name   = "Association"
  ) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  scale_edge_alpha(range = c(0.4, 0.9), guide = "none") +
  geom_node_point(
    aes(size = degree, fill = betweenness),
    shape = 21, color = "white", stroke = 0.6
  ) +
  scale_size(range = c(3, 10), name = "Degree") +
  scale_fill_gradient(
    low  = "#FEE08B", high = "#1A9850",
    name = "Betweenness\n(normalized)"
  ) +
  geom_node_text(
    aes(label = name), repel = TRUE,
    size = 2.5, color = "grey20"
  ) +
  theme_graph() +
  labs(
    title    = "Microbial Co-occurrence Network",
    subtitle = paste0("|r| ≥ ", r_threshold,
                      "  |  FDR < ", fdr_threshold,
                      "  |  Spearman, CLR-transformed")
  ) +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )
 
ggsave("cooccurrence_network.pdf", p_network,
       width = 14, height = 10, dpi = 300)
cat("Saved: cooccurrence_network.pdf\n")

# Group-stratified Analysis
# Edit "Group" to match the grouping column in your metadata
group_col <- "Group"
 
if (group_col %in% colnames(metadata)) {
 
  groups <- unique(metadata[[group_col]])
 
  cor_by_group <- lapply(groups, function(grp) {
    samps <- metadata$SampleID[metadata[[group_col]] == grp]
    samps <- intersect(samps, rownames(abund_clr))
    rcorr(as.matrix(abund_clr[samps, ]), type = "spearman")$r
  })
  names(cor_by_group) <- groups
 
  # Differential correlation (first two groups)
  if (length(groups) >= 2) {
    diff_cor <- cor_by_group[[2]] - cor_by_group[[1]]
 
    pdf("cooccurrence_differential.pdf", width = 12, height = 11)
    corrplot(
      diff_cor,
      method = "color",
      type   = "upper",
      order  = "hclust",
      col    = colorRampPalette(c("#762A83", "white", "#1B7837"))(200),
      tl.cex = 0.6,
      tl.col = "black",
      title  = paste0("Differential Co-occurrence: ",
                       groups[2], " − ", groups[1]),
      mar    = c(0, 0, 2, 0)
    )
    dev.off()
    cat("Saved: cooccurrence_differential.pdf\n")
  }
 
} else {
  cat("Column '", group_col, "' not found in metadata — skipping group analysis.\n",
      sep = "")
  cat("Available columns:", paste(colnames(metadata), collapse = ", "), "\n")
}

# Save Results
write.csv(edges,      "cooccurrence_edges.csv",      row.names = FALSE)
write.csv(node_stats, "cooccurrence_node_stats.csv", row.names = FALSE)
 
net_stats <- data.frame(
  n_nodes         = vcount(g),
  n_edges         = ecount(g),
  positive_edges  = sum(E(g)$sign == "positive"),
  negative_edges  = sum(E(g)$sign == "negative"),
  density         = round(edge_density(g), 4),
  avg_degree      = round(mean(degree(g)), 2),
  clustering_coef = round(transitivity(g, type = "global"), 4),
  top_hub         = node_stats$Taxon[1]
)

write.csv(net_stats, "cooccurrence_network_stats.csv", row.names = FALSE)
cat("\nSaved: cooccurrence_edges.csv, cooccurrence_node_stats.csv,",
    "cooccurrence_network_stats.csv\n")
