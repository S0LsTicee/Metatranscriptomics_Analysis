# Install Packages
my_packages <- c("tidyverse", "propr", "igraph", "ggraph",
              "ggrepel", "RColorBrewer", "patchwork",
              "corrplot", "scales")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(propr)
library(igraph)
library(ggraph)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(corrplot)
library(scales)

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

# Prevalence filter: keep taxa in >= 20% of samples
prev      <- rowMeans(abund_mat > 0)
taxa_keep <- names(prev[prev >= 0.20])
cat("Taxa kept after prevalence filter:", length(taxa_keep), "/",
    nrow(abund_mat), "\n")
 
abund_filt <- abund_mat[taxa_keep, ]

# Replace zeros with small pseudocount (required for propr)
abund_filt[abund_filt == 0] <- 0.5 * min(abund_filt[abund_filt > 0])

# Transpose to samples x taxa (propr expects this orientation)
count_mat <- t(abund_filt)   # samples x taxa

# Proportionality Metrics
# rho  (ρ): symmetric, -1 to 1; ρ ≈ 1 = proportional
# phi  (φ): asymmetric, 0 = perfectly proportional
# phs  (φs): symmetric version of phi
pr_rho <- propr(count_mat, metric = "rho",  ivar = "clr", symmetrize = TRUE)
pr_phi <- propr(count_mat, metric = "phs",  ivar = "clr", symmetrize = TRUE)
# Extract matrices
rho_mat <- pr_rho@matrix
phi_mat <- pr_phi@matrix
 
rownames(rho_mat) <- colnames(rho_mat) <- colnames(count_mat)
rownames(phi_mat) <- colnames(phi_mat) <- colnames(count_mat)

# FDR via permutation (updateCutoffs)
pr_rho <- updateCutoffs(pr_rho,
                         cutoffs = seq(0.5, 0.99, 0.01),
                         ncores  = 1)

# Choose rho cutoff at FDR < 0.05
fdr_table  <- pr_rho@fdr
rho_cutoff <- fdr_table$cutoff[which(fdr_table$FDR <= 0.05)[1]]
if (is.na(rho_cutoff)) {
  rho_cutoff <- 0.60
  cat("No FDR < 0.05 cutoff found — defaulting to rho >=", rho_cutoff, "\n")
} else {
  cat("Rho cutoff at FDR < 0.05:", rho_cutoff, "\n")
}

# Extract Significant Pairs
sig_pairs <- subset(pr_rho, cutoff = rho_cutoff)@pairs

if (nrow(sig_pairs) == 0) {
  cat("No pairs pass cutoff — lowering to rho >= 0.50 for illustration.\n")
  rho_cutoff <- 0.50
  sig_pairs  <- subset(pr_rho, cutoff = rho_cutoff)@pairs
}
 
# Add taxon names to pairs
sig_pairs$Taxon1 <- colnames(count_mat)[sig_pairs$Partner]
sig_pairs$Taxon2 <- colnames(count_mat)[sig_pairs$Pair]

# Per-Taxon Summary Stats
# Mean rho with all other taxa (connectivity)
rho_diag     <- rho_mat
diag(rho_diag) <- NA
taxon_summary <- data.frame(
  Taxon      = colnames(rho_mat),
  Mean_rho   = rowMeans(rho_diag, na.rm = TRUE),
  Max_rho    = apply(rho_diag, 1, max,   na.rm = TRUE),
  Min_rho    = apply(rho_diag, 1, min,   na.rm = TRUE),
  N_sig_pairs = sapply(colnames(rho_mat), function(t) {
    sum(sig_pairs$Taxon1 == t | sig_pairs$Taxon2 == t)
  }),
  stringsAsFactors = FALSE
) %>% arrange(desc(N_sig_pairs))
 
cat("\nTop taxa by number of significant proportional pairs:\n")
print(head(taxon_summary, 10))

# Rho Matrix Heatmap
pdf("proportionality_heatmap_rho.pdf", width = 12, height = 11)
corrplot(
  rho_mat,
  method   = "color",
  type     = "upper",
  order    = "hclust",
  addrect  = 4,
  col      = colorRampPalette(c("#2166AC", "white", "#D6604D"))(200),
  tl.cex   = 0.6,
  tl.col   = "black",
  title    = paste0("Proportionality Matrix (ρ)  |  cutoff = ", rho_cutoff),
  mar      = c(0, 0, 2, 0),
  cl.lim   = c(-1, 1)
)
dev.off()
cat("Saved: proportionality_heatmap_rho.pdf\n")

# Phi Matrix Heatmap
pdf("proportionality_heatmap_phi.pdf", width = 12, height = 11)
corrplot(
  phi_mat,
  method   = "color",
  type     = "upper",
  order    = "hclust",
  addrect  = 4,
  col      = colorRampPalette(c("white", "#762A83"))(200),
  tl.cex   = 0.6,
  tl.col   = "black",
  title    = "Proportionality Matrix (φs)  |  lower = more proportional",
  mar      = c(0, 0, 2, 0)
)
dev.off()
cat("Saved: proportionality_heatmap_phi.pdf\n")

# FDR Curve
fdr_df <- as.data.frame(pr_rho@fdr)
 
p_fdr <- ggplot(fdr_df, aes(x = cutoff, y = FDR)) +
  geom_line(color = "#4393C3", linewidth = 1) +
  geom_point(size = 2.5, color = "#4393C3") +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  geom_vline(xintercept = rho_cutoff, linetype = "dashed",
             color = "grey50") +
  annotate("text", x = rho_cutoff + 0.01, y = 0.5,
           label = paste0("ρ = ", rho_cutoff),
           hjust = 0, color = "grey30", size = 3.5) +
  annotate("text", x = min(fdr_df$cutoff), y = 0.06,
           label = "FDR = 0.05", hjust = 0,
           color = "#D73027", size = 3.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "FDR vs Rho Cutoff (Permutation-based)",
    subtitle = "Pairs above selected cutoff are considered significant",
    x        = "Rho (ρ) Cutoff",
    y        = "False Discovery Rate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Rho Distribution
rho_vals <- rho_mat[upper.tri(rho_mat)]
 
p_rho_dist <- ggplot(data.frame(rho = rho_vals), aes(x = rho)) +
  geom_histogram(bins = 60, fill = "#74ADD1",
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept =  rho_cutoff, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  geom_vline(xintercept = -rho_cutoff, linetype = "dashed",
             color = "#4575B4", linewidth = 0.8) +
  annotate("text", x = rho_cutoff + 0.01, y = Inf,
           label = paste0("ρ = +", rho_cutoff),
           vjust = 2, hjust = 0, color = "#D73027", size = 3) +
  annotate("text", x = -rho_cutoff - 0.01, y = Inf,
           label = paste0("ρ = -", rho_cutoff),
           vjust = 2, hjust = 1, color = "#4575B4", size = 3) +
  labs(
    title    = "Distribution of Pairwise Rho Values",
    subtitle = paste0(sum(abs(rho_vals) >= rho_cutoff),
                      " of ", length(rho_vals),
                      " pairs pass |ρ| ≥ ", rho_cutoff),
    x        = "Rho (ρ)",
    y        = "Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Top Taxa Connectivity Bar
top_conn <- taxon_summary %>%
  filter(N_sig_pairs > 0) %>%
  arrange(desc(N_sig_pairs)) %>%
  slice_head(n = 20) %>%
  mutate(Taxon = fct_reorder(Taxon, N_sig_pairs))
 
p_conn <- ggplot(top_conn,
                 aes(x = N_sig_pairs, y = Taxon, fill = Mean_rho)) +
  geom_col(alpha = 0.9, width = 0.7) +
  scale_fill_gradient2(low  = "#2166AC", mid = "grey90",
                       high = "#D6604D", midpoint = 0,
                       name = "Mean ρ") +
  labs(
    title    = paste0("Top 20 Taxa by Significant Pairs (ρ ≥ ", rho_cutoff, ")"),
    x        = "Number of Significant Proportional Pairs",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    axis.text.y        = element_text(face = "italic", size = 9),
    panel.grid.major.y = element_blank()
  )

# Proportionality Network
if (nrow(sig_pairs) > 0) {
  edges_net <- data.frame(
    from   = sig_pairs$Taxon1,
    to     = sig_pairs$Taxon2,
    rho    = sig_pairs$propr,
    stringsAsFactors = FALSE
  )
 
  g <- graph_from_data_frame(edges_net, directed = FALSE,
                              vertices = colnames(count_mat))
 
  E(g)$sign       <- ifelse(E(g)$rho > 0, "positive", "negative")
  E(g)$abs_rho    <- abs(E(g)$rho)
  V(g)$degree     <- degree(g)
  V(g)$betweenness <- betweenness(g, normalized = TRUE)
 
  # Remove isolated nodes for cleaner plot
  g_sub <- delete_vertices(g, V(g)[degree(g) == 0])
 
  set.seed(42)
  p_net <- ggraph(g_sub, layout = "fr") +
    geom_edge_link(
      aes(color = sign, width = abs_rho, alpha = abs_rho)
    ) +
    scale_edge_color_manual(
      values = c("positive" = "#D6604D", "negative" = "#2166AC"),
      name   = "Direction"
    ) +
    scale_edge_width(range = c(0.3, 2.5), guide = "none") +
    scale_edge_alpha(range = c(0.4, 0.9), guide = "none") +
    geom_node_point(
      aes(size = degree, fill = betweenness),
      shape = 21, color = "white", stroke = 0.5
    ) +
    scale_size(range = c(2, 9), name = "Degree") +
    scale_fill_gradient(low = "#FEE08B", high = "#1A9850",
                        name = "Betweenness\n(normalized)") +
    geom_node_text(
      aes(label = name), repel = TRUE,
      size = 2.4, color = "grey20", fontface = "italic"
    ) +
    theme_graph() +
    labs(
      title    = "Proportionality Network (ρ)",
      subtitle = paste0("ρ ≥ ", rho_cutoff,
                        "  |  ", ecount(g_sub), " edges  |  ",
                        vcount(g_sub), " taxa")
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey50")
    )
 
  ggsave("proportionality_network.pdf", p_net,
         width = 14, height = 10, dpi = 300)
  cat("Saved: proportionality_network.pdf\n")
}

# Two-groups group-stratified Rho 
if (n_groups == 2) {
  cat("\nRunning group-stratified proportionality...\n")
 
  pr_list <- lapply(groups, function(grp) {
    samps <- rownames(metadata[metadata[[group_col]] == grp, ])
    sub   <- count_mat[samps, ]
    pr    <- propr(sub, metric = "rho", ivar = "clr", symmetrize = TRUE)
    pr@matrix
  })
  names(pr_list) <- groups
 
  diff_rho <- pr_list[[2]] - pr_list[[1]]
 
  pdf("proportionality_differential.pdf", width = 12, height = 11)
  corrplot(
    diff_rho,
    method  = "color",
    type    = "upper",
    order   = "hclust",
    col     = colorRampPalette(c("#762A83", "white", "#1B7837"))(200),
    tl.cex  = 0.6,
    tl.col  = "black",
    title   = paste0("Differential ρ: ", groups[2], " − ", groups[1]),
    mar     = c(0, 0, 2, 0),
    cl.lim  = c(-1, 1)
  )
  dev.off()
  cat("Saved: proportionality_differential.pdf\n")
}

# Save Plots
pdf("proportionality_diagnostics.pdf", width = 14, height = 6)
print(p_fdr + p_rho_dist)
dev.off()
cat("Saved: proportionality_diagnostics.pdf\n")
 
pdf("proportionality_connectivity.pdf", width = 10, height = 7)
print(p_conn)
dev.off()
cat("Saved: proportionality_connectivity.pdf\n")

# Save Tables
write.csv(sig_pairs,     "proportionality_sig_pairs.csv",  row.names = FALSE)
write.csv(taxon_summary, "proportionality_taxon_stats.csv", row.names = FALSE)
write.csv(as.data.frame(pr_rho@fdr),
                         "proportionality_fdr_table.csv",  row.names = FALSE)
 
# Full rho matrix
write.csv(as.data.frame(rho_mat),
          "proportionality_rho_matrix.csv", row.names = TRUE)
 
cat("\nSaved tables:\n")
