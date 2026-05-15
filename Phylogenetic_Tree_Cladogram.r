# The Cladogram represents the phylogenetic tree of taxonomy-based detected species and visualizes the result in both circular and linear graphs
# Install Packages
my_packages <- c("tidyverse", "ggtree", "treeio", "ape",
                   "RColorBrewer", "patchwork", "ggtreeExtra",
                   "scales", "ggplot2")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(RColorBrewer)
library(patchwork)
library(ggtreeExtra)
library(scales)
library(ggplot2)

# Load Files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is placeholder for your merged metaphlan file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE, 
)

metadata <- read.delim(
  "metadata_merged.txt",  #This is placeholder for your merged metadata file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE
)
meta_id_col   <- colnames(metadata)[1]

# colnames(metaphlan)[1] <- "taxonomy"
# Keep species-level rows only
species_rows <- metaphlan %>%
  filter(grepl("s__", taxonomy) & !grepl("t__", taxonomy))
 
# Mean relative abundance across samples
sample_cols  <- setdiff(colnames(species_rows), "taxonomy")
species_rows <- species_rows %>%
  mutate(mean_abund = rowMeans(across(all_of(sample_cols))))
 
# Parse MetaPhlAn taxonomy string into levels
tax_parsed <- species_rows %>%
  select(taxonomy, mean_abund) %>%
  separate(taxonomy, into = c("kingdom","phylum","class","order","family","genus","species"),
           sep = "\\|", fill = "right", extra = "drop") %>%
  mutate(across(kingdom:species, ~ gsub(".__", "", .x)))  # strip k__, p__, etc. 

# Build Newick Tree from Taxonomy
build_newick_from_taxonomy <- function(tax_df) {
  # Construct nested list from taxonomy levels, then convert to phylo
  levels_used <- c("phylum","class","order","family","genus","species")
 
  # Build edge table: parent → child
  edges <- list()
  for (lvl_i in seq_along(levels_used)) {
    if (lvl_i == 1) {
      parents <- rep("root", nrow(tax_df))
    } else {
      parents <- tax_df[[levels_used[lvl_i - 1]]]
    }
    children <- tax_df[[levels_used[lvl_i]]]
    valid    <- !is.na(children) & children != ""
    edges[[lvl_i]] <- data.frame(parent = parents[valid],
                                 child  = children[valid],
                                 stringsAsFactors = FALSE)
  }
 
  edge_df <- bind_rows(edges) %>% distinct()
 
  # Unique nodes
  all_nodes   <- unique(c(edge_df$parent, edge_df$child))
  n_nodes     <- length(all_nodes)
  node_index  <- setNames(seq_along(all_nodes), all_nodes)
 
  # Tips = nodes that never appear as parent
  tips        <- setdiff(edge_df$child, edge_df$parent)
  internals   <- setdiff(all_nodes, tips)
 
  # Reorder: tips first, then internals (ape convention)
  ordered_nodes <- c(tips, internals)
  node_index    <- setNames(seq_along(ordered_nodes), ordered_nodes)
  n_tip         <- length(tips)
 
  edge_mat <- matrix(
    c(node_index[edge_df$parent], node_index[edge_df$child]),
    ncol = 2
  )
 
  phy <- list(
    edge        = edge_mat,
    tip.label   = tips,
    node.label  = internals,
    Nnode       = length(internals),
    edge.length = rep(1, nrow(edge_mat))
  )
  class(phy) <- "phylo"
  phy <- ape::reorder.phylo(phy)
  return(phy)
}
 
phy <- build_newick_from_taxonomy(tax_parsed)
message(sprintf("Tree built: %d tips, %d internal nodes.", length(phy$tip.label), phy$Nnode))

# Annotation Data Frames
# Tip annotations: phylum, mean abundance
tip_annot <- tax_parsed %>%
  select(species, phylum, mean_abund) %>%
  filter(species %in% phy$tip.label) %>%
  distinct(species, .keep_all = TRUE) %>%
  rename(label = species)
 
# Colour palette per phylum
phyla       <- sort(unique(tip_annot$phylum))
n_phyla     <- length(phyla)
phylum_pal  <- setNames(
  colorRampPalette(brewer.pal(min(n_phyla, 12), "Paired"))(n_phyla),
  phyla
)
 
tip_annot <- tip_annot %>%
  mutate(phylum_col = phylum_pal[phylum],
         log_abund  = log10(mean_abund + 1e-5))
 
# Auto-detect metadata group variable
group_col <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) & n_distinct(.) > 1 & n_distinct(.) < 20)) %>%
  colnames() %>% first()
 
# Per-group mean abundance for heatmap strip (optional)
if (!is.null(group_col) && !is.na(group_col)) {
  groups      <- unique(metadata[[group_col]])
  group_abund <- lapply(groups, function(grp) {
    samps <- metadata %>% filter(.data[[group_col]] == grp) %>% pull(meta_id_col)
    samps <- intersect(samps, sample_cols)
    if (length(samps) == 0) return(NULL)
    species_rows %>%
      select(taxonomy, all_of(samps)) %>%
      mutate(species   = gsub(".*s__", "", taxonomy),
             grp_mean  = rowMeans(across(all_of(samps)))) %>%
      select(species, grp_mean) %>%
      rename(!!grp := grp_mean)
  })
  group_abund <- group_abund[!sapply(group_abund, is.null)]
  if (length(group_abund) > 0) {
    heatmap_df <- reduce(group_abund, full_join, by = "species") %>%
      filter(species %in% phy$tip.label) %>%
      rename(label = species)
  } else {
    heatmap_df <- NULL
  }
} else {
  heatmap_df <- NULL
}

# Annotate Tree Object
annotate_tree <- function(p, tip_annot, phylum_pal, add_heatmap = FALSE) {
  p <- p %<+% tip_annot +
    geom_tippoint(aes(color = phylum, size = log_abund), alpha = 0.85) +
    scale_color_manual(values = phylum_pal, name = "Phylum", na.translate = FALSE) +
    scale_size_continuous(name = "log₁₀(Mean Abundance)", range = c(1, 5)) +
    geom_tiplab(aes(label = gsub("s__", "", label)),
                size = 2.2, offset = 0.05, align = FALSE) +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 9, face = "bold"),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9)
    )
 
  if (add_heatmap && !is.null(heatmap_df)) {
    group_cols_heat <- setdiff(colnames(heatmap_df), "label")
    heat_long <- heatmap_df %>%
      pivot_longer(-label, names_to = "group", values_to = "abundance")
    p <- p +
      new_scale_fill() +
      geom_fruit(
        data     = heat_long,
        geom     = geom_tile,
        mapping  = aes(y = label, x = group, fill = log10(abundance + 1e-5)),
        offset   = 0.08, pwidth = 0.6
      ) +
      scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#D6604D",
                           midpoint = 0, name = "log₁₀(Abund)\nby group",
                           na.value = "grey90")
  }
  return(p)
}

# Circular/Fan Cladogram
p_circ <- ggtree(phy, layout = "circular", branch.length = "none",
                  color = "grey50", linewidth = 0.35) %>%
  annotate_tree(tip_annot, phylum_pal, add_heatmap = !is.null(heatmap_df)) +
  labs(
    title    = "Taxonomy-Based Cladogram  (Circular)",
    subtitle = sprintf("%d species  |  %d phyla  |  branch length = none (cladogram)",
                       length(phy$tip.label), n_phyla)
  )
 
ggsave("cladogram_circular.pdf", p_circ,
       width = 16, height = 16, dpi = 300, limitsize = FALSE)
message("Saved: cladogram_circular.pdf")

# Rectangular/Linear Cladogram
p_rect <- ggtree(phy, layout = "rectangular", branch.length = "none",
                  color = "grey50", linewidth = 0.35) %>%
  annotate_tree(tip_annot, phylum_pal, add_heatmap = FALSE) +
  labs(
    title    = "Taxonomy-Based Cladogram  (Rectangular)",
    subtitle = sprintf("%d species  |  %d phyla", length(phy$tip.label), n_phyla)
  )
 
# Dynamically set height based on number of tips
tip_height <- max(8, length(phy$tip.label) * 0.18)
ggsave("cladogram_rectangular.pdf", p_rect,
       width = 14, height = tip_height, dpi = 300, limitsize = FALSE)
message("Saved: cladogram_rectangular.pdf")

# Phylum Legend Panel
legend_df <- data.frame(
  phylum     = names(phylum_pal),
  color      = phylum_pal,
  n_species  = as.integer(table(tip_annot$phylum)[names(phylum_pal)])
) %>% arrange(desc(n_species))
 
p_legend <- ggplot(legend_df, aes(x = reorder(phylum, n_species), y = n_species, fill = phylum)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  scale_fill_manual(values = phylum_pal) +
  coord_flip() +
  labs(title = "Species per Phylum", x = NULL, y = "Number of species") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
 
ggsave("cladogram_phylum_summary.pdf", p_legend,
       width = 8, height = max(4, n_phyla * 0.5), dpi = 300)
message("Saved: cladogram_phylum_summary.pdf")

# Save Annotation Table
write.csv(tip_annot %>% select(label, phylum, mean_abund),
          "cladogram_tip_annotations.csv", row.names = FALSE)
message("Saved: cladogram_tip_annotations.csv")
