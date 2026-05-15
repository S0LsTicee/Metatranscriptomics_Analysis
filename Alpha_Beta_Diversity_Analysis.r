# Install Packages
my_packages <- c("phyloseq", "vegan", "ggplot2", "dplyr",
  "tidyr", "readr", "RColorBrewer",
  "ggpubr", "pairwiseAdonis")

install.packages(my_packages)

# Load merged files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for the name of your merged metaphlan file
  header           = TRUE,
  sep              = "\t",
  comment.char     = "#",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.delim(
  "metadata_merged.txt", #This is a placeholder for the name of your merged metaphlan file
  header           = TRUE,
  sep              = "\t",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

# Parse MetaPhlAn Table
# Detect taxonomy column
tax_col <- grep("clade_name|#clade_name", colnames(metaphlan), value = TRUE)[1]
if (is.na(tax_col)) tax_col <- colnames(metaphlan)[1]
 
# Filter species-level rows only
species_rows  <- grepl("s__", metaphlan[[tax_col]])
metaphlan_species  <- metaphlan[species_rows, ]
message(paste("Species-level rows:", nrow(metaphlan_species)))
 
# Build numeric abundance matrix (rows = taxa, cols = samples)
taxa_labels <- metaphlan_species[[tax_col]]
abund_mat   <- as.matrix(metaphlan_species[, !colnames(metaphlan_species) %in% tax_col])
mode(abund_mat) <- "numeric"
rownames(abund_mat) <- taxa_labels

# Parse Taxonomy
parse_taxonomy <- function(clade_string) {
  ranks   <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  prefixes <- c("k__","p__","c__","o__","f__","g__","s__")
  parts   <- strsplit(clade_string, "\\|")[[1]]
  result  <- setNames(rep(NA_character_, length(ranks)), ranks)
  for (i in seq_along(prefixes)) {
    hit <- grep(paste0("^", prefixes[i]), parts, value = TRUE)
    if (length(hit)) result[ranks[i]] <- sub(prefixes[i], "", hit[1])
  }
  result
}
 
tax_df <- do.call(rbind, lapply(taxa_labels, parse_taxonomy))
rownames(tax_df) <- taxa_labels

# Align Samples 
sample_id_col  <- colnames(metadata)[1]
rownames(metadata) <- metadata[[sample_id_col]]
 
common_samples <- intersect(colnames(abund_mat), rownames(metadata))
message(paste("Common samples:", length(common_samples)))
if (length(common_samples) == 0)
  stop("No matching sample IDs. Check column names vs metadata first column.")
 
abund_mat <- abund_mat[, common_samples, drop = FALSE]
metadata  <- metadata[common_samples, , drop = FALSE]
 
# Detect primary grouping variable (first non-ID column)
group_var <- colnames(metadata)[colnames(metadata) != sample_id_col][1]
message(paste("Grouping variable:", group_var))

# Build Phyloseq Object
OTU <- otu_table(abund_mat, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_df))
META <- sample_data(metadata)
ps  <- phyloseq(OTU, TAX, META)
message(ps)

# Alpha Diversity Plot 
# Calculate metrics
alpha_df <- estimate_richness(ps, measures = c(
  "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"
))
alpha_df$SampleID <- rownames(alpha_df)

# Pielou's Evenness (Shannon / log(Observed))
alpha_df$Evenness <- alpha_df$Shannon / log(alpha_df$Observed)

# Merge with metadata
alpha_merged <- merge(alpha_df, data.frame(sample_data(ps)),
                      by.x = "SampleID", by.y = sample_id_col)
 
write.csv(alpha_merged, "alpha_diversity_results.csv", row.names = FALSE)
message("Saved: alpha_diversity_results.csv")

# Box Plot
metrics_plot <- c("Observed","Chao1","Shannon","Simpson","InvSimpson","Evenness")
metrics_plot <- intersect(metrics_plot, colnames(alpha_merged))
 
alpha_long <- pivot_longer(
  alpha_merged,
  cols      = all_of(metrics_plot),
  names_to  = "Metric",
  values_to = "Value"
)
 
p_alpha <- ggplot(alpha_long,
                  aes(x = .data[[group_var]], y = Value,
                      fill = .data[[group_var]])) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2,
               alpha = 0.75, width = 0.55) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, shape = 16) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title  = "Alpha Diversity",
       x      = group_var,
       y      = "Value",
       fill   = group_var) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92", color = "grey60"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )
 
ggsave("alpha_diversity_boxplots.pdf", p_alpha, width = 14, height = 10)
message("Saved: alpha_diversity_boxplots.pdf")

# Statistical tests — Kruskal-Wallis + Wilcoxon pairwise
kw_list <- lapply(metrics_plot, function(m) {
  formula <- as.formula(paste(m, "~", group_var))
  kw      <- kruskal.test(formula, data = alpha_merged)
  data.frame(Metric    = m,
             KW_stat   = round(kw$statistic, 4),
             df        = kw$parameter,
             KW_pvalue = round(kw$p.value,   6))
})
kw_df <- do.call(rbind, kw_list)
write.csv(kw_df, "alpha_kruskal_wallis.csv", row.names = FALSE)
message("Saved: alpha_kruskal_wallis.csv")
print(kw_df)

# Beta Diversity Matrice
# Proportional abundance matrix — samples x taxa
relab_mat <- t(otu_table(ps)) / 100   # convert MetaPhlAn % → proportions

# Distance matrices
dist_methods <- c("bray", "jaccard", "euclidean")
 
dist_list <- lapply(setNames(dist_methods, dist_methods), function(m) {
  if (m == "jaccard") {
    vegdist(relab_mat, method = "jaccard", binary = TRUE)
  } else {
    vegdist(relab_mat, method = m)
  }
})

# Aitchison (CLR-based Euclidean) — robust for compositional data
clr_transform <- function(mat) {
  # mat: samples x taxa (proportions)
  mat[mat == 0] <- 1e-6          # pseudocount
  log_mat <- log(mat)
  sweep(log_mat, 1, rowMeans(log_mat), "-")
}
clr_mat           <- clr_transform(as.matrix(relab_mat))
dist_list$aitchison <- dist(clr_mat, method = "euclidean")
 
all_methods <- names(dist_list)

# Save distance matrices
for (m in all_methods) {
  write.csv(as.matrix(dist_list[[m]]),
            paste0("beta_", m, "_matrix.csv"))
}
message("Beta distance matrices saved.")

# PCoA Ordination
run_pcoa <- function(dist_obj, label, group_col, meta_df, outfile) {
  pcoa   <- cmdscale(dist_obj, eig = TRUE, k = 2)
  eig    <- pcoa$eig
  var_ex <- round(eig / sum(eig[eig > 0]) * 100, 1)
 
  pcoa_df <- data.frame(
    PC1      = pcoa$points[, 1],
    PC2      = pcoa$points[, 2],
    SampleID = rownames(pcoa$points)
  )
  pcoa_df <- merge(pcoa_df, meta_df,
                   by.x = "SampleID", by.y = sample_id_col)
 
  p <- ggplot(pcoa_df,
              aes(x = PC1, y = PC2, color = .data[[group_col]])) +
    geom_point(size = 4, alpha = 0.85) +
    stat_ellipse(aes(group = .data[[group_col]]),
                 level = 0.95, linetype = 2, linewidth = 0.7) +
    scale_color_brewer(palette = "Set1") +
    labs(title  = paste("PCoA –", label),
         x      = paste0("PC1 (", var_ex[1], "%)"),
         y      = paste0("PC2 (", var_ex[2], "%)"),
         color  = group_col) +
    theme_bw(base_size = 13) +
    theme(
      plot.title     = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
 
  ggsave(outfile, p, width = 8, height = 6)
  message(paste("Saved:", outfile))
  invisible(pcoa_df)
}
 
run_pcoa(dist_list$bray,      "Bray-Curtis",       group_var, metadata, "pcoa_braycurtis.pdf")
run_pcoa(dist_list$jaccard,   "Jaccard (binary)",  group_var, metadata, "pcoa_jaccard.pdf")
run_pcoa(dist_list$aitchison, "Aitchison (CLR)",   group_var, metadata, "pcoa_aitchison.pdf")
run_pcoa(dist_list$euclidean, "Euclidean",         group_var, metadata, "pcoa_euclidean.pdf")

# NMDA (Bray-Curtis)
set.seed(42)
nmds <- metaMDS(relab_mat, distance = "bray", k = 2, trymax = 100, trace = FALSE)
message(paste("NMDS stress:", round(nmds$stress, 4)))
 
nmds_df <- data.frame(
  NMDS1    = nmds$points[, 1],
  NMDS2    = nmds$points[, 2],
  SampleID = rownames(nmds$points)
)
nmds_df <- merge(nmds_df, metadata,
                 by.x = "SampleID", by.y = sample_id_col)
 
p_nmds <- ggplot(nmds_df,
                 aes(x = NMDS1, y = NMDS2, color = .data[[group_var]])) +
  geom_point(size = 4, alpha = 0.85) +
  stat_ellipse(aes(group = .data[[group_var]]),
               level = 0.95, linetype = 2, linewidth = 0.7) +
  scale_color_brewer(palette = "Set1") +
  annotate("text", x = Inf, y = -Inf,
           label = paste("Stress:", round(nmds$stress, 4)),
           hjust = 1.1, vjust = -0.5, size = 4, color = "grey40") +
  labs(title = "NMDS – Bray-Curtis",
       color  = group_var) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
 
ggsave("nmds_braycurtis.pdf", p_nmds, width = 8, height = 6)
message("Saved: nmds_braycurtis.pdf")

# PERMANOVA
permanova_results <- lapply(all_methods, function(m) {
  formula <- as.formula(paste("dist_list[[m]] ~", group_var))
  result  <- adonis2(formula, data = data.frame(metadata), permutations = 999)
  cat("\nPERMANOVA –", m, "\n")
  print(result)
  data.frame(
    Distance  = m,
    R2        = round(result$R2[1],         4),
    F_stat    = round(result$F[1],          4),
    p_value   = round(result$`Pr(>F)`[1],   4)
  )
})
 
perm_df <- do.call(rbind, permanova_results)
write.csv(perm_df, "permanova_results.csv", row.names = FALSE)
message("Saved: permanova_results.csv")
print(perm_df)

# PERMDISP
groups <- factor(metadata[[group_var]])
 
permdisp_results <- lapply(all_methods, function(m) {
  disp <- betadisper(dist_list[[m]], groups)
  test <- permutest(disp, permutations = 999)
  cat("\nPERMDISP –", m, "\n")
  print(test$tab)
  data.frame(
    Distance = m,
    F_stat   = round(test$tab$F[1], 4),
    p_value  = round(test$tab$`Pr(>F)`[1], 4)
  )
})
 
permdisp_df <- do.call(rbind, permdisp_results)
write.csv(permdisp_df, "permdisp_results.csv", row.names = FALSE)
message("Saved: permdisp_results.csv")
print(permdisp_df)

# Pairwise PERMANOVA
if (requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  message("\n--- Pairwise PERMANOVA (Bray-Curtis, BH correction) ---")
  pw <- pairwiseAdonis::pairwise.adonis(
    dist_list$bray,
    factors    = groups,
    p.adjust.m = "BH",
    perm       = 999
  )
  print(pw)
  write.csv(pw, "pairwise_permanova_braycurtis.csv", row.names = FALSE)
  message("Saved: pairwise_permanova_braycurtis.csv")
}
