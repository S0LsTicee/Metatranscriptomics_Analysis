#Install Packages
my_packages <- c("tidyverse", "vegan", "ggplot2", "ggrepel", "RColorBrewer", "patchwork")

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

#Load Files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE,
)

metadata <- read.delim(
  "metadata_merged.txt", #This is a placeholder for your merged metaphlan file
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE
)

# MetaPhlAn clade names: keep rows that contain "s__" but NOT "t__"
species_data <- metaphlan %>%
  filter(grepl("s__", .[[1]]) & !grepl("t__", .[[1]]))

# Rename taxonomy column for clarity
colnames(species_data)[1] <- "taxonomy"
 
# Extract short species name (last element after "|")
species_data <- species_data %>%
  mutate(species = gsub(".*\\|", "", taxonomy)) %>%
  select(species, everything(), -taxonomy)
 
# Convert to long format
species_long <- species_data %>%
  pivot_longer(-species, names_to = "sample_id", values_to = "abundance")

# Identify the sample-ID column in metadata (first column by default)
meta_id_col <- colnames(metadata)[1]
 
species_long <- species_long %>%
  left_join(metadata, by = setNames(meta_id_col, "sample_id"))

# Rank Abundance Per Sample
rank_abund <- species_long %>%
  filter(abundance > 0) %>%
  group_by(sample_id) %>%
  arrange(desc(abundance), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Mean Rank-Abundance Across All Samples
mean_rank_abund <- rank_abund %>%
  group_by(rank) %>%
  summarise(
    mean_abundance = mean(abundance),
    sd_abundance   = sd(abundance),
    n              = n(),
    se             = sd_abundance / sqrt(n),
    .groups = "drop"
  )

# Most abundant species labels (top 10 by mean)
top_species <- rank_abund %>%
  group_by(species) %>%
  summarise(mean_ab = mean(abundance), .groups = "drop") %>%
  slice_max(mean_ab, n = 10) %>%
  pull(species)
 
label_df <- rank_abund %>%
  filter(species %in% top_species) %>%
  group_by(species) %>%
  summarise(
    rank           = median(rank),
    mean_abundance = mean(abundance),
    .groups = "drop"
  )

# Mean Rank-Abundance Curve
p_mean <- ggplot(mean_rank_abund, aes(x = rank, y = log10(mean_abundance + 1e-5))) +
  geom_ribbon(aes(ymin = log10(pmax(mean_abundance - se, 1e-5)),
                  ymax = log10(mean_abundance + se + 1e-5)),
              fill = "#4E9AF1", alpha = 0.25) +
  geom_line(color = "#1565C0", linewidth = 0.9) +
  geom_point(color = "#1565C0", size = 1.5, alpha = 0.7) +
  geom_text_repel(
    data        = label_df,
    aes(x = rank, y = log10(mean_abundance + 1e-5), label = gsub("s__", "", species)),
    size        = 3,
    max.overlaps = 15,
    box.padding = 0.4
  ) +
  labs(
    title    = "Mean Rank-Abundance Curve",
    subtitle = "Shaded band = ±1 SE  |  Top 10 species labelled",
    x        = "Rank",
    y        = "log₁₀(Mean Relative Abundance)"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# Per-group rank-abundance (if group col exists)
# Auto-detect a usable grouping variable (non-numeric, >1 but <20 levels)
group_col <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) && n_distinct(.) > 1 && n_distinct(.) < 20)) %>%
  colnames() %>%
  first()
 
if (!is.null(group_col) && !is.na(group_col)) {
  message(sprintf("Grouping by column: '%s'", group_col))
 
  group_rank <- rank_abund %>%
    group_by(.data[[group_col]], rank) %>%
    summarise(mean_abundance = mean(abundance), .groups = "drop")
 
  n_groups  <- n_distinct(group_rank[[group_col]])
  pal        <- colorRampPalette(brewer.pal(min(n_groups, 8), "Dark2"))(n_groups)
 
  p_group <- ggplot(group_rank,
                    aes(x = rank,
                        y = log10(mean_abundance + 1e-5),
                        color = .data[[group_col]])) +
    geom_line(linewidth = 0.85, alpha = 0.85) +
    scale_color_manual(values = pal) +
    labs(
      title    = paste("Rank-Abundance by", group_col),
      x        = "Rank",
      y        = "log₁₀(Mean Relative Abundance)",
      color    = group_col
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold"),
      legend.position = "right"
    )
 
  combined_plot <- p_mean / p_group
  ggsave("rank_abundance_plots.pdf", combined_plot,
         width = 10, height = 12, dpi = 300)
  message("Saved: rank_abundance_plots.pdf")
 
} else {
  message("No suitable grouping column detected; saving mean-curve plot only.")
  ggsave("rank_abundance_plots.pdf", p_mean,
         width = 10, height = 6, dpi = 300)
  message("Saved: rank_abundance_plots.pdf")
}

# Per-sample rank-abundance facet plot (up to 20 samples)
samples_to_plot <- unique(rank_abund$sample_id)
if (length(samples_to_plot) > 20) {
  message("More than 20 samples – plotting first 20 in facet grid.")
  samples_to_plot <- samples_to_plot[1:20]
}
 
p_facet <- rank_abund %>%
  filter(sample_id %in% samples_to_plot) %>%
  ggplot(aes(x = rank, y = log10(abundance + 1e-5))) +
  geom_line(color = "#E53935", linewidth = 0.6, alpha = 0.8) +
  geom_point(size = 0.8, alpha = 0.6, color = "#B71C1C") +
  facet_wrap(~ sample_id, scales = "free_x") +
  labs(
    title = "Per-Sample Rank-Abundance Curves",
    x     = "Rank",
    y     = "log₁₀(Relative Abundance)"
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title   = element_text(face = "bold"),
    strip.text   = element_text(size = 7),
    axis.text    = element_text(size = 7)
  )
 
ggsave("rank_abundance_per_sample.pdf", p_facet,
       width = 14, height = 10, dpi = 300)
message("Saved: rank_abundance_per_sample.pdf")

# Save Ranked Tables
write.csv(rank_abund,       "rank_abundance_long.csv",       row.names = FALSE)
write.csv(mean_rank_abund,  "rank_abundance_mean_summary.csv", row.names = FALSE)
message("Saved: rank_abundance_long.csv")
message("Saved: rank_abundance_mean_summary.csv")
 
message("\n✓ Rank-abundance analysis complete.")
