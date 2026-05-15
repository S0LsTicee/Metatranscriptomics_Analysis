# Install Packages
my_packages <- c("tidyverse", "ggplot2", "ggrepel", "RColorBrewer",
              "patchwork", "scales", "phyloseq", "vegan")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(scales)
library(phyloseq)
library(vegan)

#Load Files
metaphlan <- read.table(
  "metaphlan_merged.txt", #This is the placeholder for your merged metaphlan file
  header       = TRUE,
  sep          = "\t",
  row.names    = 1,
  check.names  = FALSE,
)

metadata <- read.table(
  "metadata_merged.txt", #This is the placeholder for your merged metaphlan file
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

# Overall prevalence-abundance statistics
n_samples <- ncol(abund_mat)
 
taxon_stats <- data.frame(
  Taxon        = rownames(abund_mat),
  Prevalence   = rowMeans(abund_mat > 0) * 100,        # %
  MeanAbund    = rowMeans(abund_mat),                   # %
  MedianAbund  = apply(abund_mat, 1, median),
  SDAbund      = apply(abund_mat, 1, sd),
  MaxAbund     = apply(abund_mat, 1, max),
  N_present    = rowSums(abund_mat > 0),
  stringsAsFactors = FALSE
) %>%
  mutate(
    CV          = SDAbund / (MeanAbund + 1e-9) * 100,   # coefficient of variation
    Category    = case_when(
      Prevalence >= 80 & MeanAbund >= 1   ~ "High prev / High abund",
      Prevalence >= 80 & MeanAbund <  1   ~ "High prev / Low abund",
      Prevalence <  80 & MeanAbund >= 1   ~ "Low prev / High abund",
      TRUE                                ~ "Low prev / Low abund"
    )
  ) %>%
  arrange(desc(Prevalence), desc(MeanAbund))
 
cat("\nTaxon category summary:\n")
print(table(taxon_stats$Category))

# Per-group Prevalence Abundance
group_stats <- lapply(groups, function(grp) {
  samps   <- rownames(metadata[metadata[[group_col]] == grp, ])
  sub_mat <- abund_mat[, samps, drop = FALSE]
  data.frame(
    Taxon       = rownames(sub_mat),
    Group       = grp,
    Prevalence  = rowMeans(sub_mat > 0) * 100,
    MeanAbund   = rowMeans(sub_mat),
    MedianAbund = apply(sub_mat, 1, median),
    SDAbund     = apply(sub_mat, 1, sd),
    N_present   = rowSums(sub_mat > 0),
    N_total     = ncol(sub_mat),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Differential Prevalence Between Two Groups 
if (n_groups == 2) {
  prev_wide <- group_stats %>%
    select(Taxon, Group, Prevalence) %>%
    pivot_wider(names_from = Group, values_from = Prevalence,
                names_prefix = "Prev_")
 
  abund_wide <- group_stats %>%
    select(Taxon, Group, MeanAbund) %>%
    pivot_wider(names_from = Group, values_from = MeanAbund,
                names_prefix = "Abund_")
 
  diff_df <- left_join(prev_wide, abund_wide, by = "Taxon") %>%
    mutate(
      Prev_diff  = .data[[paste0("Prev_",  groups[2])]] -
                   .data[[paste0("Prev_",  groups[1])]],
      Abund_diff = .data[[paste0("Abund_", groups[2])]] -
                   .data[[paste0("Abund_", groups[1])]]
    ) %>%
    arrange(desc(abs(Prev_diff)))
 
  cat("\nTop taxa by differential prevalence (", groups[2], " - ",
      groups[1], "):\n", sep = "")
  print(head(diff_df %>% select(Taxon, Prev_diff, Abund_diff), 10))
}

# Overall prevalence vs abundance plot (4-quadrant) 
prev_cutoff  <- 50    # %
abund_cutoff <- 0.1   # %
 
label_taxa <- taxon_stats %>%
  filter(Prevalence >= prev_cutoff | MeanAbund >= 1) %>%
  arrange(desc(Prevalence)) %>%
  slice_head(n = 25)
 
p1 <- ggplot(taxon_stats,
             aes(x = MeanAbund + 1e-4, y = Prevalence,
                 color = Category, size = MeanAbund)) +
  geom_vline(xintercept = abund_cutoff, linetype = "dashed",
             color = "grey60", linewidth = 0.6) +
  geom_hline(yintercept = prev_cutoff,  linetype = "dashed",
             color = "grey60", linewidth = 0.6) +
  geom_point(alpha = 0.75) +
  geom_text_repel(
    data = label_taxa,
    aes(label = Taxon), size = 2.4, color = "grey20",
    fontface = "italic", max.overlaps = 25, show.legend = FALSE
  ) +
  scale_x_log10(
    labels = label_log(digits = 2),
    breaks = 10^seq(-4, 2, 1)
  ) +
  scale_color_manual(
    values = c(
      "High prev / High abund" = "#D73027",
      "High prev / Low abund"  = "#FC8D59",
      "Low prev / High abund"  = "#4575B4",
      "Low prev / Low abund"   = "#ABD9E9"
    ),
    name = "Category"
  ) +
  scale_size(range = c(1, 6), guide = "none") +
  annotate("text", x = 0.001, y = 95, label = "Rare biosphere",
           color = "grey50", size = 3, hjust = 0) +
  annotate("text", x = 5,    y = 95, label = "Core microbiome",
           color = "#D73027", size = 3, fontface = "bold", hjust = 0) +
  labs(
    title    = "Prevalence vs Mean Relative Abundance",
    subtitle = paste0("n = ", n_samples, " samples  |  ",
                      nrow(taxon_stats), " species"),
    x        = "Mean Relative Abundance (%, log scale)",
    y        = "Prevalence (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey50"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# Prevalence distribution histogram
p2 <- ggplot(taxon_stats, aes(x = Prevalence)) +
  geom_histogram(binwidth = 5, fill = "#4393C3",
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = prev_cutoff, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  annotate("text", x = prev_cutoff + 1, y = Inf,
           label = paste0(prev_cutoff, "% cutoff"),
           vjust = 2, hjust = 0, color = "#D73027", size = 3.5) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(
    title = "Prevalence Distribution",
    x     = "Prevalence (%)",
    y     = "Number of Taxa"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Abundance distribution histogram
p3 <- ggplot(taxon_stats %>% filter(MeanAbund > 0),
             aes(x = MeanAbund)) +
  geom_histogram(bins = 40, fill = "#74ADD1",
                 color = "white", alpha = 0.85) +
  scale_x_log10(labels = label_log(digits = 2)) +
  geom_vline(xintercept = abund_cutoff, linetype = "dashed",
             color = "#D73027", linewidth = 0.8) +
  annotate("text", x = abund_cutoff * 1.5, y = Inf,
           label = paste0(abund_cutoff, "% cutoff"),
           vjust = 2, hjust = 0, color = "#D73027", size = 3.5) +
  labs(
    title = "Mean Abundance Distribution",
    x     = "Mean Relative Abundance (%, log scale)",
    y     = "Number of Taxa"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Top 20 taxa by prevalence bar graph
top20_prev <- taxon_stats %>%
  arrange(desc(Prevalence)) %>%
  slice_head(n = 20) %>%
  mutate(Taxon = fct_reorder(Taxon, Prevalence))
 
p4 <- ggplot(top20_prev,
             aes(x = Prevalence, y = Taxon, fill = MeanAbund)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_vline(xintercept = prev_cutoff, linetype = "dashed",
             color = "grey50") +
  scale_fill_gradient(low = "#DEEBF7", high = "#08519C",
                      name = "Mean\nAbundance (%)") +
  scale_x_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  labs(
    title = "Top 20 Taxa by Prevalence",
    x     = "Prevalence (%)",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    axis.text.y        = element_text(face = "italic", size = 9),
    panel.grid.major.y = element_blank()
  )

# Top 20 taxa by mean abundance bar graph
top20_abund <- taxon_stats %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 20) %>%
  mutate(Taxon = fct_reorder(Taxon, MeanAbund))
 
p5 <- ggplot(top20_abund,
             aes(x = MeanAbund, y = Taxon, fill = Prevalence)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_errorbarh(
    aes(xmin = pmax(MeanAbund - SDAbund, 0),
        xmax = MeanAbund + SDAbund),
    height = 0.3, color = "grey40", linewidth = 0.4
  ) +
  scale_fill_gradient(low = "#FEE5D9", high = "#A50F15",
                      name = "Prevalence (%)") +
  labs(
    title    = "Top 20 Taxa by Mean Abundance",
    subtitle = "Error bars = ± 1 SD",
    x        = "Mean Relative Abundance (%)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold"),
    axis.text.y        = element_text(face = "italic", size = 9),
    panel.grid.major.y = element_blank()
  )

# Per-group prevalence vs abundance faceted
top_taxa_labels <- taxon_stats %>%
  filter(Prevalence >= prev_cutoff) %>%
  pull(Taxon)
 
p6 <- ggplot(group_stats,
             aes(x = MeanAbund + 1e-4, y = Prevalence,
                 color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    data = group_stats %>% filter(Taxon %in% top_taxa_labels),
    aes(label = Taxon), size = 2.2, fontface = "italic",
    max.overlaps = 15, show.legend = FALSE
  ) +
  geom_hline(yintercept = prev_cutoff, linetype = "dashed",
             color = "grey60") +
  scale_x_log10(labels = label_log(digits = 2)) +
  scale_color_manual(values = pal, name = group_col) +
  facet_wrap(~ Group, ncol = 2) +
  labs(
    title = "Prevalence vs Abundance by Group",
    x     = "Mean Relative Abundance (%, log scale)",
    y     = "Prevalence (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "none"
  )

# Differential prevalence plot between two groups
if (n_groups == 2 && exists("diff_df")) {
  diff_plot_df <- diff_df %>%
    arrange(desc(abs(Prev_diff))) %>%
    slice_head(n = 25) %>%
    mutate(Taxon     = fct_reorder(Taxon, Prev_diff),
           Direction = ifelse(Prev_diff > 0, groups[2], groups[1]))
 
  p7 <- ggplot(diff_plot_df,
               aes(x = Prev_diff, y = Taxon, fill = Direction)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_vline(xintercept = 0, color = "grey40") +
    scale_fill_manual(values = pal[1:2], name = "Higher in") +
    labs(
      title    = paste0("Differential Prevalence: ",
                         groups[2], " vs ", groups[1]),
      subtitle = "Top 25 taxa by absolute prevalence difference",
      x        = paste0("Prevalence difference (%, ",
                         groups[2], " − ", groups[1], ")"),
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold"),
      axis.text.y        = element_text(face = "italic", size = 9),
      panel.grid.major.y = element_blank()
    )
}

# Save Plots
pdf("prev_abund_main.pdf", width = 12, height = 8)
print(p1)
dev.off()
cat("Saved: prev_abund_main.pdf\n")
 
pdf("prev_abund_distributions.pdf", width = 12, height = 5)
print(p2 + p3)
dev.off()
cat("Saved: prev_abund_distributions.pdf\n")
 
pdf("prev_abund_top_taxa.pdf", width = 14, height = 6)
print(p4 + p5)
dev.off()
cat("Saved: prev_abund_top_taxa.pdf\n")
 
pdf("prev_abund_by_group.pdf",
    width = 6 * min(n_groups, 2), height = 5 * ceiling(n_groups / 2))
print(p6)
dev.off()
cat("Saved: prev_abund_by_group.pdf\n")
 
if (exists("p7")) {
  pdf("prev_abund_differential.pdf", width = 10, height = 8)
  print(p7)
  dev.off()
  cat("Saved: prev_abund_differential.pdf\n")
}

# Save Tables
write.csv(taxon_stats,  "prev_abund_overall.csv",   row.names = FALSE)
write.csv(group_stats,  "prev_abund_by_group.csv",  row.names = FALSE)
 
if (exists("diff_df")) {
  write.csv(diff_df, "prev_abund_differential.csv", row.names = FALSE)
  cat("Saved: prev_abund_differential.csv\n")
}
