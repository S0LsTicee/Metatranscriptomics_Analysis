# Install packages 
my_packages <- c("phyloseq", "vegan", "ggplot2", "dplyr",
  "tidyr", "RColorBrewer", "iNEXT")

# Check for missing files and install them
not_installed <- not_installed[!(not_installed %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

# Load files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header           = TRUE,
  sep              = "\t",
  comment.char     = "#",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)
 
metadata <- read.delim(
  "metadata_merged.txt", #This is a placeholder for your merged metadata file
  header           = TRUE,
  sep              = "\t",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

# Parse metaphlan table
tax_col <- grep("clade_name|#clade_name", colnames(metaphlan), value = TRUE)[1]
if (is.na(tax_col)) tax_col <- colnames(metaphlan)[1]
 
# Species-level only
species_rows <- grepl("s__", metaphlan[[tax_col]])
metaphlan_species <- metaphlan[species_rows, ]
message(paste("Species-level rows:", nrow(metaphlan_species)))
 
taxa_labels <- metaphlan_species[[tax_col]]
abund_mat   <- as.matrix(metaphlan_species[, !colnames(metaphlan_species) %in% tax_col])
mode(abund_mat) <- "numeric"
rownames(abund_mat) <- taxa_labels

# Align Samples
sample_id_col       <- colnames(metadata)[1]
rownames(metadata)  <- metadata[[sample_id_col]]
common_samples      <- intersect(colnames(abund_mat), rownames(metadata))
 
if (length(common_samples) == 0)
  stop("No matching sample IDs. Check sample names in both files.")
 
abund_mat <- abund_mat[, common_samples, drop = FALSE]
metadata  <- metadata[common_samples, , drop = FALSE]
 
group_var <- colnames(metadata)[colnames(metadata) != sample_id_col][1]

# Convert Abundances to Integer Counts
# MetaPhlAn outputs relative abundances (%).
# For rarefaction, we scale to integer pseudo-counts (sum = 10,000 per sample).
scale_depth <- 10000
count_mat   <- round(abund_mat / 100 * scale_depth)   # taxa x samples
count_mat   <- count_mat[rowSums(count_mat) > 0, ]    # remove all-zero taxa
 
# Sample read depths (column sums)
sample_depths <- colSums(count_mat)
message("Sample depths (pseudo-counts):")
print(summary(sample_depths))

# Build Phyloseq Object
OTU  <- otu_table(count_mat, taxa_are_rows = TRUE)
META <- sample_data(metadata)
ps   <- phyloseq(OTU, META)
message(ps)

# Choose Rarefaction Depth
# Standard: minimum sample depth (drops no samples).
# Adjust rarefy_depth manually if you want to drop low-depth outliers.
 
rarefy_depth <- min(sample_depths)
message(paste("Rarefaction depth (min sample depth):", rarefy_depth))
 
# Warn if any sample is far below the median
median_depth <- median(sample_depths)
low_samples  <- names(sample_depths[sample_depths < median_depth * 0.5])
if (length(low_samples) > 0) {
  message("WARNING: Low-depth samples detected (< 50% of median):")
  print(low_samples)
  message("Consider raising rarefy_depth and excluding these samples.")
}

# Rarefaction Curves (iNEXT)
# iNEXT expects a list or matrix: species x samples, integer counts
inext_input <- as.data.frame(count_mat)   # taxa x samples
 
set.seed(42)
inext_out <- iNEXT(
  inext_input,
  q         = 0,        # Hill number q=0 → species richness
  datatype  = "abundance",
  knots     = 40,
  se        = TRUE,
  conf      = 0.95,
  nboot     = 100
)
 
# Extract rarefaction data
rare_df <- do.call(rbind, lapply(names(inext_out$iNextEst), function(s) {
  df        <- inext_out$iNextEst[[s]]
  df$Sample <- s
  df
}))
 
# Merge group info
rare_meta <- data.frame(
  Sample = rownames(metadata),
  Group  = metadata[[group_var]],
  stringsAsFactors = FALSE
)
rare_df <- merge(rare_df, rare_meta, by = "Sample")
 
# 8a. Rarefaction curve plot — coloured by sample
n_samples <- length(common_samples)
palette   <- if (n_samples <= 12) {
  colorRampPalette(brewer.pal(min(n_samples, 12), "Paired"))(n_samples)
} else {
  colorRampPalette(brewer.pal(8, "Set2"))(n_samples)
}
 
p_rare_sample <- ggplot(rare_df,
                        aes(x = m, y = qD, group = Sample, color = Sample)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Sample),
              alpha = 0.12, color = NA) +
  geom_vline(xintercept = rarefy_depth,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text",
           x     = rarefy_depth * 1.02,
           y     = max(rare_df$qD.UCL, na.rm = TRUE) * 0.95,
           label = paste("Rarefaction\ndepth =", rarefy_depth),
           hjust = 0, size = 3.5, color = "red") +
  scale_color_manual(values = palette) +
  scale_fill_manual(values  = palette) +
  labs(title    = "Rarefaction Curves – Per Sample",
       subtitle = paste("q = 0 (Species Richness) | Depth =", scale_depth, "pseudo-counts"),
       x        = "Sampling Depth (pseudo-counts)",
       y        = "Species Richness",
       color    = "Sample", fill = "Sample") +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = if (n_samples > 20) "none" else "right"
  )
 
ggsave("rarefaction_curves_per_sample.pdf", p_rare_sample,
       width = 12, height = 7)
message("Saved: rarefaction_curves_per_sample.pdf")

# Rarefaction curve plot — coloured by group (mean ± CI per group)
group_rare <- rare_df %>%
  group_by(Group, m) %>%
  summarise(
    mean_rich = mean(qD),
    se_rich   = sd(qD) / sqrt(n()),
    .groups   = "drop"
  )
 
n_groups   <- length(unique(group_rare$Group))
grp_palette <- brewer.pal(max(3, min(n_groups, 8)), "Set1")[seq_len(n_groups)]
 
p_rare_group <- ggplot(group_rare,
                       aes(x = m, y = mean_rich,
                           color = Group, fill = Group)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = mean_rich - se_rich,
                  ymax = mean_rich + se_rich),
              alpha = 0.2, color = NA) +
  geom_vline(xintercept = rarefy_depth,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  annotate("text",
           x     = rarefy_depth * 1.02,
           y     = max(group_rare$mean_rich + group_rare$se_rich,
                       na.rm = TRUE) * 0.95,
           label = paste("Rarefaction\ndepth =", rarefy_depth),
           hjust = 0, size = 3.5, color = "red") +
  scale_color_manual(values = grp_palette) +
  scale_fill_manual(values  = grp_palette) +
  labs(title    = paste("Rarefaction Curves – By", group_var),
       subtitle = "Mean ± SE across samples within group",
       x        = "Sampling Depth (pseudo-counts)",
       y        = "Species Richness",
       color    = group_var, fill = group_var) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "right"
  )
 
ggsave("rarefaction_curves_by_group.pdf", p_rare_group,
       width = 10, height = 6)
message("Saved: rarefaction_curves_by_group.pdf")

# Rarefy the OTU Table
  ps_rare <- rarefy_even_depth(
  ps,
  sample.size = rarefy_depth,
  rngseed     = 42,
  replace     = FALSE,
  trimOTUs    = TRUE,
  verbose     = TRUE
)
 
message(ps_rare)
 
# Save rarefied count table
rare_count_table <- as.data.frame(otu_table(ps_rare))
write.csv(rare_count_table, "rarefied_count_table.csv")
message("Saved: rarefied_count_table.csv")

# Alpha Diversity Rarefied 
alpha_rare <- estimate_richness(
ps_rare,
measures = c("Observed","Shannon","Simpson","InvSimpson","Chao1")
)
alpha_rare$Evenness <- alpha_rare$Shannon / log(alpha_rare$Observed)
alpha_rare$SampleID <- rownames(alpha_rare)
 
alpha_rare_merged <- merge(alpha_rare, data.frame(sample_data(ps_rare)),
                           by.x = "SampleID", by.y = sample_id_col)
 
write.csv(alpha_rare_merged, "alpha_diversity_rarefied.csv", row.names = FALSE)
message("Saved: alpha_diversity_rarefied.csv")

# # Alpha Diversity Rarefied Box Plots
metrics_plot <- c("Observed","Shannon","Simpson","InvSimpson","Evenness")
metrics_plot <- intersect(metrics_plot, colnames(alpha_rare_merged))
 
alpha_long <- pivot_longer(alpha_rare_merged,
                           cols      = all_of(metrics_plot),
                           names_to  = "Metric",
                           values_to = "Value")
 
p_alpha <- ggplot(alpha_long,
                  aes(x = .data[[group_var]], y = Value,
                      fill = .data[[group_var]])) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2,
               alpha = 0.75, width = 0.55) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.7) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title    = "Alpha Diversity (Rarefied)",
       subtitle = paste("Rarefaction depth:", rarefy_depth),
       x        = group_var, y = "Value", fill = group_var) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    legend.position  = "bottom"
  )
 
ggsave("alpha_diversity_rarefied_boxplots.pdf", p_alpha,
       width = 14, height = 10)
message("Saved: alpha_diversity_rarefied_boxplots.pdf")
 
# Kruskal-Wallis on rarefied alpha
kw_list <- lapply(metrics_plot, function(m) {
  kw <- kruskal.test(as.formula(paste(m, "~", group_var)),
                     data = alpha_rare_merged)
  data.frame(Metric   = m,
             KW_stat  = round(kw$statistic, 4),
             df       = kw$parameter,
             p_value  = round(kw$p.value, 6))
})
kw_df <- do.call(rbind, kw_list)
write.csv(kw_df, "alpha_kruskal_wallis_rarefied.csv", row.names = FALSE)
message("Saved: alpha_kruskal_wallis_rarefied.csv")
print(kw_df)

# Depth Summary
depth_df <- data.frame(
  SampleID = names(sample_depths),
  Depth    = sample_depths,
  Group    = metadata[names(sample_depths), group_var]
)
 
p_depth <- ggplot(depth_df,
                  aes(x = reorder(SampleID, -Depth),
                      y = Depth, fill = Group)) +
  geom_bar(stat = "identity", color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = rarefy_depth,
             color = "red", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = nrow(depth_df) * 0.6,
           y = rarefy_depth * 1.04,
           label = paste("Rarefaction depth =", rarefy_depth),
           color = "red", size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Sample Sequencing Depth",
       x     = "Sample", y = "Pseudo-count Depth",
       fill  = group_var) +
  theme_bw(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )
 
ggsave("sample_depth_barplot.pdf", p_depth,
       width = max(10, n_samples * 0.35), height = 6)
message("Saved: sample_depth_barplot.pdf")
 
write.csv(depth_df, "sample_depth_summary.csv", row.names = FALSE)
message("Saved: sample_depth_summary.csv")
