#This code uses DESeq2 and ALDEx2 to visualize fold change and significance on the taxonomy tree
#Install Packages
my_packages <- c("tidyverse", "ggtree", "ape", "DESeq2",
                   "ALDEx2", "RColorBrewer", "scales",
                   "ggplot2", "patchwork", "ggrepel")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(ggtree)
library(ape)
library(DESeq2)
library(ALDEx2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(patchwork)
library(ggrepel)

#Load Files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE, 
)

metadata <- read.delim(
  "metadata_merged.txt",  #This is a placeholder for your merged metadata file
  header = TRUE,
  sep = "\t", check.names = FALSE
)
meta_id_col   <- colnames(metadata)[1]

group_col <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) & n_distinct(.) > 1 & n_distinct(.) < 20)) %>%
  colnames() %>% first()
 
if (is.null(group_col) || is.na(group_col))
  stop("No suitable grouping variable found. Set group_col manually.")
message(sprintf("Grouping variable: '%s'  |  Levels: %s",
                group_col, paste(unique(metadata[[group_col]]), collapse = ", ")))

# Build Count Matrix colnames(metaphlan)[1] <- "taxonomy"
species_rows <- metaphlan %>%
  filter(grepl("s__", taxonomy) & !grepl("t__", taxonomy)) %>%
  mutate(species = gsub(".*s__", "s__", taxonomy))
 
sample_cols <- intersect(colnames(species_rows), metadata[[meta_id_col]])
 
abund_mat <- species_rows %>%
  select(species, all_of(sample_cols)) %>%
  column_to_rownames("species") %>%
  as.data.frame()
 
# Convert relative abundance → pseudo-counts (multiply by 1e6, round)
count_mat <- round(abund_mat * 1e6)
count_mat <- count_mat[rowSums(count_mat) > 0, ]
 
# Align metadata to sample order
meta_aligned <- metadata[match(sample_cols, metadata[[meta_id_col]]), ]
meta_aligned[[group_col]] <- factor(meta_aligned[[group_col]])

# DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta_aligned,
  design    = as.formula(paste("~", group_col))
)
dds <- DESeq(dds, quiet = TRUE)
 
# Extract results for all pairwise contrasts
groups    <- levels(meta_aligned[[group_col]])
contrasts <- combn(groups, 2, simplify = FALSE)
 
deseq_all <- lapply(contrasts, function(pr) {
  res <- tryCatch(
    results(dds, contrast = c(group_col, pr[1], pr[2]),
            alpha = 0.05, independentFiltering = TRUE),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  as.data.frame(res) %>%
    rownames_to_column("species") %>%
    mutate(contrast = paste(pr[1], "vs", pr[2]),
           ref = pr[2], comparison = pr[1])
}) %>% bind_rows()
 
deseq_sig <- deseq_all %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, comparison, ref))

# ALDEx2 Validation & Pairwise
aldex_results <- tryCatch({
  conds <- as.character(meta_aligned[[group_col]])
  if (length(unique(conds)) == 2) {
    ax <- ALDEx2::aldex(count_mat, conds, mc.samples = 128,
                        test = "t", effect = TRUE, verbose = FALSE)
    ax %>% rownames_to_column("species") %>%
      select(species, we.ep, we.eBH, effect) %>%
      rename(aldex_p = we.ep, aldex_padj = we.eBH, aldex_effect = effect)
  } else {
    message("ALDEx2 skipped for >2 groups (use Kruskal mode manually).")
    NULL
  }
}, error = function(e) { message("ALDEx2 error: ", e$message); NULL })
 
# Join DESeq2 + ALDEx2
if (!is.null(aldex_results)) {
  deseq_sig <- deseq_sig %>%
    left_join(aldex_results, by = "species") %>%
    mutate(validated = !is.na(aldex_padj) & aldex_padj < 0.05)
} else {
  deseq_sig$validated <- NA
}

# Build Taxonomy Tree
tax_parsed <- species_rows %>%
  select(taxonomy, species) %>%
  separate(taxonomy,
           into = c("kingdom","phylum","class","order","family","genus","sp2"),
           sep = "\\|", fill = "right", extra = "drop") %>%
  mutate(across(everything(), ~ gsub(".__", "", .x)),
         species = sp2) %>% select(-sp2) %>%
  distinct(species, .keep_all = TRUE)
 
build_tree <- function(tax_df) {
  levels_used <- c("phylum","class","order","family","genus","species")
  edges <- list()
  for (i in seq_along(levels_used)) {
    parents  <- if (i == 1) rep("root", nrow(tax_df)) else tax_df[[levels_used[i-1]]]
    children <- tax_df[[levels_used[i]]]
    valid    <- !is.na(children) & children != ""
    edges[[i]] <- data.frame(parent = parents[valid], child = children[valid])
  }
  edge_df   <- bind_rows(edges) %>% distinct()
  tips      <- setdiff(edge_df$child, edge_df$parent)
  internals <- setdiff(unique(c(edge_df$parent, edge_df$child)), tips)
  all_nodes <- c(tips, internals)
  node_idx  <- setNames(seq_along(all_nodes), all_nodes)
  edge_mat  <- matrix(c(node_idx[edge_df$parent], node_idx[edge_df$child]), ncol = 2)
  phy <- list(edge = edge_mat, tip.label = tips, node.label = internals,
              Nnode = length(internals), edge.length = rep(1, nrow(edge_mat)))
  class(phy) <- "phylo"
  ape::reorder.phylo(phy)
}
 
phy <- build_tree(tax_parsed)

# Tip Annotation
# Pick best contrast per species (lowest padj)
best_hit <- deseq_sig %>%
  group_by(species) %>%
  slice_min(padj, n = 1) %>%
  ungroup() %>%
  select(species, log2FoldChange, padj, direction, contrast, validated)
 
dir_groups  <- unique(best_hit$direction)
dir_pal     <- setNames(
  colorRampPalette(brewer.pal(min(length(dir_groups), 8), "Set1"))(length(dir_groups)),
  dir_groups
)
 
tip_annot <- tax_parsed %>%
  rename(label = species) %>%
  left_join(best_hit, by = c("label" = "species")) %>%
  mutate(
    is_sig   = !is.na(padj),
    pt_size  = ifelse(is_sig, rescale(abs(log2FoldChange), to = c(2, 7)), 1.2),
    pt_color = ifelse(is_sig, direction, NA_character_)
  )

# Circular Cladogram
p_circ <- ggtree(phy, layout = "circular", branch.length = "none",
                  color = "grey55", linewidth = 0.3) %<+% tip_annot +
  geom_tippoint(aes(color = pt_color, size = pt_size), alpha = 0.88) +
  geom_tiplab(aes(label = ifelse(is_sig, gsub("s__","", label), "")),
              size = 1.9, offset = 0.08) +
  scale_color_manual(values = dir_pal, na.value = "grey80",
                     name = "Enriched in", na.translate = FALSE) +
  scale_size_identity() +
  labs(
    title    = "Differential Abundance Cladogram  (Circular)",
    subtitle = sprintf("DESeq2 padj < 0.05  |  %d significant species  |  Colour = enriched group",
                       n_distinct(deseq_sig$species))
  ) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9),
        legend.position = "right")
 
ggsave("cladogram_diffabund_circular.pdf", p_circ,
       width = 16, height = 16, dpi = 300, limitsize = FALSE)
message("Saved: cladogram_diffabund_circular.pdf")

# Rectangular Cladogram
p_rect <- ggtree(phy, layout = "rectangular", branch.length = "none",
                  color = "grey55", linewidth = 0.3) %<+% tip_annot +
  geom_tippoint(aes(color = pt_color, size = pt_size), alpha = 0.88) +
  geom_tiplab(aes(label = gsub("s__","", label)), size = 2.1, offset = 0.05) +
  scale_color_manual(values = dir_pal, na.value = "grey80",
                     name = "Enriched in", na.translate = FALSE) +
  scale_size_identity() +
  labs(title    = "Differential Abundance Cladogram  (Rectangular)",
       subtitle = sprintf("%d significant species", n_distinct(deseq_sig$species))) +
  theme_tree2() +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")
 
tip_h <- max(8, length(phy$tip.label) * 0.18)
ggsave("cladogram_diffabund_rectangular.pdf", p_rect,
       width = 14, height = tip_h, dpi = 300, limitsize = FALSE)
message("Saved: cladogram_diffabund_rectangular.pdf")

# Volcano Plot
p_volc <- deseq_all %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig   = padj < 0.05 & abs(log2FoldChange) > 1,
    label = ifelse(sig, gsub("s__","", species), NA_character_),
    color = case_when(
      sig & log2FoldChange > 0 ~ comparison,
      sig & log2FoldChange < 0 ~ ref,
      TRUE                     ~ "ns"
    )
  ) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = color, label = label)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(size = 2.5, max.overlaps = 20, box.padding = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", color = "grey40") +
  scale_color_manual(
    values = c(setNames(dir_pal, names(dir_pal)), ns = "grey70"),
    name   = "Enriched in"
  ) +
  facet_wrap(~ contrast, scales = "free_x") +
  labs(title = "Differential Abundance Volcano Plot",
       x = "log₂ Fold Change", y = "-log₁₀(adj. p-value)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
 
ggsave("cladogram_diffabund_volcano.pdf", p_volc,
       width = 12, height = 5 * ceiling(length(contrasts) / 2),
       dpi = 300, limitsize = FALSE)
message("Saved: cladogram_diffabund_volcano.pdf")

# Save Graphs
write.csv(deseq_all  %>% arrange(padj), "deseq2_results.csv",    row.names = FALSE)
write.csv(deseq_sig  %>% arrange(padj), "deseq2_significant.csv", row.names = FALSE)
if (!is.null(aldex_results))
  write.csv(aldex_results %>% arrange(aldex_padj), "aldex2_results.csv", row.names = FALSE)
 
message(sprintf("\n✓ Differential abundance cladogram complete.  Significant species: %d",
                n_distinct(deseq_sig$species)))
