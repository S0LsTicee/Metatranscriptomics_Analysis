# The script makes an LEfSe-style cladogram that identifies taxa enriched per group via the Kruskal-Wallis method.
# Install Packages
my_packages <- c("tidyverse", "ggtree", "ape", "RColorBrewer",
                   "MASS", "scales", "ggplot2", "patchwork")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(ggtree)
library(ape)
library(RColorBrewer)
library(MASS)
library(scales)
library(ggplot2)
library(patchwork)

# Load Files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for your merged metaphlan file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE, 
)

metadata <- read.delim(
  "metadata_merged.txt", #This is a placeholder for your merged metadata file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE
)
meta_id_col   <- colnames(metadata)[1]

# Auto-detect grouping variable
group_col <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) & n_distinct(.) > 1 & n_distinct(.) < 20)) %>%
  colnames() %>% first()
 
if (is.null(group_col) || is.na(group_col))
  stop("No suitable grouping variable found in metadata. Please set group_col manually.")

# Filter Species & Build Abundance Matrix
colnames(metaphlan_raw)[1] <- "taxonomy"
species_rows <- metaphlan %>%
  filter(grepl("s__", taxonomy) & !grepl("t__", taxonomy)) %>%
  mutate(species = gsub(".*s__", "s__", taxonomy))
 
sample_cols <- intersect(colnames(species_rows), metadata[[meta_id_col]])
abund_mat   <- species_rows %>%
  select(species, all_of(sample_cols)) %>%
  column_to_rownames("species") %>%
  as.data.frame()

#LEfSe Statistical Pipeline
# Kruskal-Wallis across all groups
kw_results <- apply(abund_mat, 1, function(x) {
  df <- data.frame(abund = as.numeric(x),
                   group = metadata[[group_col]][match(sample_cols, metadata[[meta_id_col]])])
  tryCatch(kruskal.test(abund ~ group, data = df)$p.value, error = function(e) NA)
})
 
kw_df <- data.frame(species = rownames(abund_mat), kw_p = kw_results) %>%
  filter(!is.na(kw_p)) %>%
  mutate(kw_p_adj = p.adjust(kw_p, method = "BH"))
 
sig_species <- kw_df %>% filter(kw_p_adj < 0.05) %>% pull(species)
message(sprintf("Kruskal-Wallis significant species: %d / %d",
                length(sig_species), nrow(abund_mat)))

# Pairwise Wilcoxon to assign enriched group
groups      <- unique(metadata[[group_col]])
group_pairs <- combn(groups, 2, simplify = FALSE)
 
wilcox_results <- lapply(sig_species, function(sp) {
  x      <- as.numeric(abund_mat[sp, ])
  grp    <- metadata[[group_col]][match(sample_cols, metadata[[meta_id_col]])]
  means  <- tapply(x, grp, mean)
  enrich <- names(which.max(means))
 
  pvals <- sapply(group_pairs, function(pr) {
    a <- x[grp == pr[1]]; b <- x[grp == pr[2]]
    if (length(a) < 2 || length(b) < 2) return(1)
    tryCatch(wilcox.test(a, b)$p.value, error = function(e) 1)
  })
 
  data.frame(species = sp, enriched_group = enrich,
             wilcox_min_p = min(pvals), stringsAsFactors = FALSE)
})
 
wilcox_df <- bind_rows(wilcox_results) %>%
  filter(wilcox_min_p < 0.05)

# LDA score (log10 of effect size via MASS::lda)
compute_lda_score <- function(sp) {
  x   <- as.numeric(abund_mat[sp, ]) + 1e-6
  grp <- metadata[[group_col]][match(sample_cols, metadata[[meta_id_col]])]
  df  <- data.frame(abund = x, group = grp)
  tryCatch({
    fit   <- lda(group ~ abund, data = df)
    score <- abs(fit$scaling[1]) * log10(mean(x) * 1000 + 1)
    score
  }, error = function(e) NA_real_)
}
 
wilcox_df$lda_score <- sapply(wilcox_df$species, compute_lda_score)
lefse_results       <- wilcox_df %>%
  filter(!is.na(lda_score) & lda_score > 2) %>%   # LDA score threshold = 2
  arrange(desc(lda_score))

# Build Taxonomy Tree
tax_parsed <- species_rows %>%
  select(taxonomy) %>%
  mutate(species = gsub(".*s__", "s__", taxonomy)) %>%
  filter(species %in% rownames(abund_mat)) %>%
  separate(taxonomy,
           into = c("kingdom","phylum","class","order","family","genus","species2"),
           sep = "\\|", fill = "right", extra = "drop") %>%
  mutate(across(everything(), ~ gsub(".__", "", .x)),
         species = species2) %>%
  select(-species2) %>%
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
  edge_df    <- bind_rows(edges) %>% distinct()
  tips       <- setdiff(edge_df$child, edge_df$parent)
  internals  <- setdiff(unique(c(edge_df$parent, edge_df$child)), tips)
  all_nodes  <- c(tips, internals)
  node_idx   <- setNames(seq_along(all_nodes), all_nodes)
  edge_mat   <- matrix(c(node_idx[edge_df$parent], node_idx[edge_df$child]), ncol = 2)
  phy <- list(edge = edge_mat, tip.label = tips, node.label = internals,
              Nnode = length(internals), edge.length = rep(1, nrow(edge_mat)))
  class(phy) <- "phylo"
  ape::reorder.phylo(phy)
}
 
phy <- build_tree(tax_parsed)

# Tip  Node Annotation
groups_pal <- setNames(
  colorRampPalette(brewer.pal(min(length(groups), 8), "Set1"))(length(groups)),
  groups
)
 
tip_annot <- tax_parsed %>%
  mutate(label = species) %>%
  left_join(lefse_results %>% select(species, enriched_group, lda_score),
            by = "species") %>%
  mutate(
    is_sig    = !is.na(enriched_group),
    node_col  = ifelse(is_sig, groups_pal[enriched_group], "grey70"),
    node_size = ifelse(is_sig, scales::rescale(lda_score, to = c(2, 7)), 1.5)
  )

# Circular LEfSe Cladogram
p_circ <- ggtree(phy, layout = "circular", branch.length = "none",
                  color = "grey60", linewidth = 0.3) %<+% tip_annot +
  geom_tippoint(aes(color = enriched_group, size = node_size), alpha = 0.9) +
  geom_tiplab(aes(label = ifelse(is_sig, gsub("s__","", label), "")),
              size = 2.0, offset = 0.1) +
  scale_color_manual(values = groups_pal, na.value = "grey80",
                     name = group_col, na.translate = FALSE) +
  scale_size_identity() +
  labs(
    title    = "LEfSe-Style Cladogram  (Circular)",
    subtitle = sprintf("LDA score > 2  |  BH-adj Kruskal-Wallis p < 0.05  |  %d enriched taxa",
                       nrow(lefse_results))
  ) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9),
        legend.position = "right")
 
ggsave("cladogram_lefse_circular.pdf", p_circ,
       width = 16, height = 16, dpi = 300, limitsize = FALSE)
message("Saved: cladogram_lefse_circular.pdf")

# Rectangular Cladogram
p_rect <- ggtree(phy, layout = "rectangular", branch.length = "none",
                  color = "grey60", linewidth = 0.3) %<+% tip_annot +
  geom_tippoint(aes(color = enriched_group, size = node_size), alpha = 0.9) +
  geom_tiplab(aes(label = gsub("s__","", label)), size = 2.2, offset = 0.05) +
  scale_color_manual(values = groups_pal, na.value = "grey80",
                     name = group_col, na.translate = FALSE) +
  scale_size_identity() +
  labs(title    = "LEfSe-Style Cladogram  (Rectangular)",
       subtitle = sprintf("%d enriched taxa  |  LDA > 2", nrow(lefse_results))) +
  theme_tree2() +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right")
 
tip_h <- max(8, length(phy$tip.label) * 0.18)
ggsave("cladogram_lefse_rectangular.pdf", p_rect,
       width = 14, height = tip_h, dpi = 300, limitsize = FALSE)
message("Saved: cladogram_lefse_rectangular.pdf")

# LDA Bar Chart
p_lda <- lefse_results %>%
  slice_max(lda_score, n = 30) %>%
  mutate(species = gsub("s__","", species),
         species = fct_reorder(species, lda_score)) %>%
  ggplot(aes(x = lda_score, y = species, fill = enriched_group)) +
  geom_col(alpha = 0.85, color = "white") +
  scale_fill_manual(values = groups_pal, name = group_col) +
  labs(title = "Top 30 LEfSe Taxa by LDA Score",
       x = "LDA Score (log10)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
 
ggsave("cladogram_lefse_lda_bar.pdf", p_lda,
       width = 10, height = 9, dpi = 300)
message("Saved: cladogram_lefse_lda_bar.pdf")

# Save Graph
write.csv(lefse_results, "lefse_results.csv", row.names = FALSE)
message("Saved: lefse_results.csv")
