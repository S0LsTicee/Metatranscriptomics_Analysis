# Install Packages
my_packages <- c("tidyverse", "microbiome", "phyloseq", "ggplot2",
              "ggrepel", "RColorBrewer", "patchwork", "vegan",
              "reshape2")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(vegan)
library(reshape2)

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
colnames(metadata)[1] <- "SampleID"
rownames(metadata)    <- metadata$SampleID

# Species-level only
metaphlan <- metaphlan[grepl("s__", rownames(metaphlan)), ]
 
# Shorten rownames to species epithet for readability
rownames(metaphlan) <- gsub(".*s__", "", rownames(metaphlan))

# Build Phyloseq Object
common_samples <- intersect(colnames(metaphlan),
                             metadata$SampleID)
cat("Common samples:", length(common_samples), "\n")
 
otu  <- otu_table(as.matrix(metaphlan[, common_samples]),
                  taxa_are_rows = TRUE)
meta <- sample_data(metadata[common_samples, ])
 
ps <- phyloseq(otu, meta)
cat("phyloseq object created:", ntaxa(ps), "taxa,",
    nsamples(ps), "samples\n")

# Set grouping variable
group_col <- "Group"
 
if (!group_col %in% colnames(metadata)) {
  cat("Warning: '", group_col, "' not found — using 'All'.\n", sep = "")
  sample_data(ps)$Group <- "All"
  group_col <- "Group"
}
 
groups     <- unique(sample_data(ps)[[group_col]])
n_groups   <- length(groups)
pal        <- if (n_groups <= 8) brewer.pal(max(3, n_groups), "Set2") else
                colorRampPalette(brewer.pal(8, "Set2"))(n_groups)
names(pal) <- groups

# Overall Core Microbiome
# prevalence thresholds to test
prev_thresholds <- c(0.50, 0.70, 0.80, 0.90, 1.00)
abund_threshold <- 0.01   # minimum mean relative abundance (%)

# Prevalence & mean abundance per taxon
abund_mat  <- as.data.frame(otu_table(ps))   # taxa x samples
n_samp     <- ncol(abund_mat)

taxon_stats <- data.frame(
Taxon      = rownames(abund_mat),
Prevalence = rowMeans(abund_mat > 0),
MeanAbund  = rowMeans(abund_mat),
MedianAbund= apply(abund_mat, 1, median),
stringsAsFactors = FALSE
) %>% arrange(desc(Prevalence), desc(MeanAbund))

# Core at 50% prevalence & 0.01% mean abundance 
core_taxa <- taxon_stats %>%
  filter(Prevalence >= 0.50, MeanAbund >= abund_threshold) %>%
  pull(Taxon)
print(core_taxa)

# Core size across prevalence thresholds
core_sizes <- sapply(prev_thresholds, function(p) {
  sum(taxon_stats$Prevalence >= p & taxon_stats$MeanAbund >= abund_threshold)
})
core_threshold_df <- data.frame(
  Prevalence_threshold = prev_thresholds * 100,
  N_core_taxa          = core_sizes
)
print(core_threshold_df)

# Group-specific core microbiomes
core_by_group <- lapply(groups, function(grp) {
  samps  <- rownames(sample_data(ps)[sample_data(ps)[[group_col]] == grp, ])
  sub_mat <- abund_mat[, samps, drop = FALSE]
  data.frame(
    Taxon      = rownames(sub_mat),
    Prevalence = rowMeans(sub_mat > 0),
    MeanAbund  = rowMeans(sub_mat),
    Group      = grp,
    stringsAsFactors = FALSE
  )
})
core_by_group_df <- bind_rows(core_by_group)

# Core per group at 50% prevalence
core_per_group <- core_by_group_df %>%
  filter(Prevalence >= 0.50, MeanAbund >= abund_threshold) %>%
  group_by(Group) %>%
  summarise(N_core = n(), .groups = "drop")

print(core_per_group)
 
# Group-specific vs shared core
core_lists <- lapply(groups, function(grp) {
  core_by_group_df %>%
    filter(Group == grp,
           Prevalence >= 0.50,
           MeanAbund  >= abund_threshold) %>%
    pull(Taxon)
})
names(core_lists) <- groups
 
shared_core <- Reduce(intersect, core_lists)
if (length(shared_core) > 0) print(shared_core)

# Overall prevalence vs abundance scatter 
p1 <- ggplot(taxon_stats,
             aes(x = Prevalence * 100, y = MeanAbund,
                 color = Prevalence >= 0.50 & MeanAbund >= abund_threshold,
                 size  = MeanAbund)) +
  geom_point(alpha = 0.75) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = abund_threshold, linetype = "dashed",
             color = "grey50") +
  geom_text_repel(
    data = taxon_stats %>%
      filter(Prevalence >= 0.50, MeanAbund >= abund_threshold),
    aes(label = Taxon), size = 2.5, color = "grey20",
    fontface = "italic", max.overlaps = 20, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("TRUE" = "#D6604D", "FALSE" = "#92C5DE"),
    labels = c("TRUE" = "Core", "FALSE" = "Non-core"),
    name   = NULL
  ) +
  scale_size(range = c(1.5, 6), guide = "none") +
  scale_y_log10() +
  labs(
    title    = "Prevalence vs Mean Relative Abundance",
    subtitle = "Core = ≥50% prevalence & ≥0.01% mean abundance",
    x        = "Prevalence (%)",
    y        = "Mean Relative Abundance (%, log scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "top"
  )

# Core size vs prevalence threshold
p2 <- ggplot(core_threshold_df,
             aes(x = Prevalence_threshold, y = N_core_taxa)) +
  geom_line(color = "#4393C3", linewidth = 1) +
  geom_point(size = 4, color = "#4393C3", fill = "white",
             shape = 21, stroke = 1.5) +
  geom_text(aes(label = N_core_taxa), vjust = -1, size = 3.5,
            color = "grey30") +
  scale_x_continuous(breaks = prev_thresholds * 100) +
  labs(
    title = "Core Microbiome Size by Prevalence Threshold",
    x     = "Prevalence Threshold (%)",
    y     = "Number of Core Taxa"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Core taxa heatmap (prevalence per group)
heatmap_df <- core_by_group_df %>%
  filter(Taxon %in% core_taxa) %>%
  mutate(Taxon = fct_reorder(Taxon, Prevalence, .fun = mean))
 
p3 <- ggplot(heatmap_df,
             aes(x = Group, y = Taxon, fill = Prevalence * 100)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", Prevalence * 100)),
            size = 3, color = "grey20") +
  scale_fill_gradient(
    low  = "#DEEBF7", high = "#08519C",
    name = "Prevalence (%)", limits = c(0, 100)
  ) +
  labs(
    title    = "Core Taxa Prevalence by Group",
    subtitle = "Core defined at ≥50% overall prevalence",
    x        = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text.y     = element_text(face = "italic", size = 9),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

# Mean abundance of core taxa (bar chart) 
core_abund_df <- taxon_stats %>%
  filter(Taxon %in% core_taxa) %>%
  mutate(Taxon = fct_reorder(Taxon, MeanAbund))
 
p4 <- ggplot(core_abund_df,
             aes(x = MeanAbund, y = Taxon, fill = Prevalence * 100)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_errorbarh(
    aes(xmin = MeanAbund - MedianAbund * 0.1,
        xmax = MeanAbund + MedianAbund * 0.1),
    height = 0.3, color = "grey40"
  ) +
  scale_fill_gradient(low = "#9ECAE1", high = "#08306B",
                      name = "Prevalence (%)") +
  labs(
    title = "Mean Relative Abundance of Core Taxa",
    x     = "Mean Relative Abundance (%)",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(face = "italic", size = 9),
    panel.grid.major.y = element_blank()
  )

#Group-stratified prevalence per core taxon
p5 <- ggplot(
    heatmap_df,
    aes(x = Prevalence * 100, y = Taxon, color = Group)
  ) +
  geom_point(size = 4, alpha = 0.85) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = pal, name = group_col) +
  labs(
    title    = "Core Taxa Prevalence per Group",
    subtitle = "Dashed line = 50% threshold",
    x        = "Prevalence (%)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(face = "italic", size = 9),
    panel.grid.minor = element_blank()
  )

# Save Plots
pdf("core_microbiome_overview.pdf", width = 14, height = 7)
print(p1 + p2)
dev.off()
cat("Saved: core_microbiome_overview.pdf\n")
 
pdf("core_microbiome_heatmap.pdf", width = 10, height = 8)
print(p3)
dev.off()
cat("Saved: core_microbiome_heatmap.pdf\n")
 
pdf("core_microbiome_abundance.pdf", width = 12, height = 6)
print(p4 + p5)
dev.off()
cat("Saved: core_microbiome_abundance.pdf\n")

# Save Tables
write.csv(taxon_stats,
          "core_taxon_stats.csv", row.names = FALSE)
 
write.csv(
  taxon_stats %>% filter(Taxon %in% core_taxa),
  "core_taxa_overall.csv", row.names = FALSE
)
 
write.csv(core_by_group_df,
          "core_taxa_by_group.csv", row.names = FALSE)
 
write.csv(
  data.frame(Taxon = shared_core),
  "core_taxa_shared.csv", row.names = FALSE
)
 
write.csv(core_threshold_df,
          "core_threshold_summary.csv", row.names = FALSE)
 
cat("\nSaved tables:\n")
