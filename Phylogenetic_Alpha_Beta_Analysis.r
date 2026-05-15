# Install Packages
my_packages <- c("phyloseq", "vegan", "ape", "picante",
  "ggplot2", "dplyr", "tidyr", "readr",
  "RColorBrewer", "pairwiseAdonis", "GUniFrac")
install.packages(my_packages)

# Check for missing files and install them
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(phyloseq)
library(vegan)
library(ape)
library(picante)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(pairwiseAdonis)
library(GUniFrac)

#Load merged files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder name for your merged metaphlan file 
  header     = TRUE,
  sep        = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Load Merged Files
metadata <- read.delim(
  "metadata_merged.txt", #This is a placeholder name for your merged metadata file
  header     = TRUE,
  sep        = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Detect the taxonomy column (usually "clade_name" or "#clade_name")
tax_col <- grep("clade_name|#clade_name", colnames(metaphlan), value = TRUE)[1]
if (is.na(tax_col)) tax_col <- colnames(metaphlan)[1]

# Filter to species-level (rows containing 's__')
species_rows <- grepl("s__", metaphlan[[tax_col]])
metaphlan_species <- metaphlan[species_rows, ]

# Build abundance matrix (numeric, samples as columns)
taxa_names_raw <- metaphlan_species[[tax_col]]
abund_mat <- as.matrix(metaphlan_species[, !colnames(metaphlan_species) %in% tax_col])
class(abund_mat) <- "numeric"
rownames(abund_mat) <- taxa_names_raw

# Parse taxonomy into Ranks
parse_taxonomy <- function(clade_string) {
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  prefixes <- c("k__","p__","c__","o__","f__","g__","s__")
  parts <- strsplit(clade_string, "\\|")[[1]]
  result <- setNames(rep(NA_character_, length(ranks)), ranks)
  for (i in seq_along(prefixes)) {
    hit <- grep(paste0("^", prefixes[i]), parts, value = TRUE)
    if (length(hit)) result[ranks[i]] <- sub(prefixes[i], "", hit[1])
  }
  result
}

tax_table_df <- do.call(rbind, lapply(taxa_names_raw, parse_taxonomy))
rownames(tax_table_df) <- taxa_names_raw

# Align Samples Between Metadata & Abundance Tables
# Identify sample ID column in metadata (first column)
sample_id_col <- colnames(metadata)[1]
rownames(metadata) <- metadata[[sample_id_col]]

common_samples <- intersect(colnames(abund_mat), rownames(metadata))
 
if (length(common_samples) == 0) {
  stop("No matching sample IDs found. Check that sample names in metaphlan_merged.txt columns match the first column of metadata_merged.txt.")
}
 
abund_mat  <- abund_mat[, common_samples, drop = FALSE]
metadata   <- metadata[common_samples, , drop = FALSE]

#Build Phyloseq Object
OTU  <- otu_table(abund_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(as.matrix(tax_table_df))
META <- sample_data(metadata)
 
ps <- phyloseq(OTU, TAX, META)
message(ps)

# Construct an approximate Phylogenetic Tree using ape::as.phylo 
#If you have a reference tree, load it with ape::read.tree("your_tree.nwk")
build_tax_tree <- function(tax_df) {
  # Create Newick-compatible labels from species column
  sp_labels <- make.unique(
    ifelse(!is.na(tax_df[, "Species"]),
           gsub("[^A-Za-z0-9_]", "_", tax_df[, "Species"]),
           paste0("sp_", seq_len(nrow(tax_df))))
  )

  # Build classification matrix for ape
  class_df <- data.frame(
    Kingdom = tax_df[, "Kingdom"],
    Phylum  = tax_df[, "Phylum"],
    Class   = tax_df[, "Class"],
    Order   = tax_df[, "Order"],
    Family  = tax_df[, "Family"],
    Genus   = tax_df[, "Genus"],
    Species = sp_labels,
    stringsAsFactors = FALSE
  )
  class_df[is.na(class_df)] <- "Unknown"
  rownames(class_df) <- sp_labels
  ape::as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus/Species,
                data = class_df, collapse = FALSE)
}
 
phy_tree <- tryCatch({
  tree <- build_tax_tree(tax_table_df)
  # Match tip labels to OTU rownames
  rownames(abund_mat) <- make.unique(
    ifelse(!is.na(tax_table_df[, "Species"]),
           gsub("[^A-Za-z0-9_]", "_", tax_table_df[, "Species"]),
           paste0("sp_", seq_len(nrow(tax_table_df))))
  )
  tree
}, error = function(e) {
  message("Tree construction failed: ", e$message)
  message("Proceeding without a phylogenetic tree (Faith's PD will be skipped).")
  NULL
})
 
# Rebuild phyloseq with updated row names and tree
OTU  <- otu_table(abund_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(as.matrix(tax_table_df))
 
if (!is.null(phy_tree)) {
  ps <- phyloseq(OTU, TAX, META, phy_tree(phy_tree))
} else {
  ps <- phyloseq(OTU, TAX, META)
}

# Alpha Diversity Plot
# Standard richness & diversity metrics
alpha_std <- estimate_richness(ps, measures = c("Observed","Shannon","Simpson","InvSimpson","Chao1"))
alpha_std$SampleID <- rownames(alpha_std)

# Faith's Phylogenetic Alpha Diversity 
if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
  comm_mat <- t(otu_table(ps))   # samples x species
  faith_pd <- pd(comm_mat, phy_tree(ps), include.root = FALSE)
  alpha_std$Faith_PD  <- faith_pd[rownames(alpha_std), "PD"]
  alpha_std$SR_picante <- faith_pd[rownames(alpha_std), "SR"]
  message("Faith's PD calculated.")
} else {
  message("Skipping Faith's PD (no tree available).")
}
 
# Merge with metadata
alpha_df <- merge(alpha_std, data.frame(sample_data(ps)),
                  by.x = "SampleID", by.y = sample_id_col)
 
write.csv(alpha_df, "alpha_diversity_results.csv", row.names = FALSE)
message("Alpha diversity saved to: alpha_diversity_results.csv")

# Detect the first grouping variable in metadata (skip SampleID column)
group_var <- colnames(metadata)[colnames(metadata) != sample_id_col][1]
 
metrics_to_plot <- intersect(
  c("Observed","Shannon","Simpson","InvSimpson","Chao1","Faith_PD"),
  colnames(alpha_df)
)
 
alpha_long <- pivot_longer(
  alpha_df,
  cols      = all_of(metrics_to_plot),
  names_to  = "Metric",
  values_to = "Value"
)
 
p_alpha <- ggplot(alpha_long, aes(x = .data[[group_var]], y = Value, fill = .data[[group_var]])) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Alpha Diversity Metrics", x = group_var, y = "Value") +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey90"),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )
 
ggsave("alpha_diversity_boxplots.pdf", p_alpha, width = 14, height = 10)
message("Alpha diversity plot saved to: alpha_diversity_boxplots.pdf")

# Kruskal-Wallis Test for Alpha Diversity
kw_results <- lapply(metrics_to_plot, function(m) {
  formula  <- as.formula(paste(m, "~", group_var))
  test     <- kruskal.test(formula, data = alpha_df)
  data.frame(Metric = m, Statistic = round(test$statistic, 4),
             df = test$parameter, p_value = round(test$p.value, 6))
})
kw_df <- do.call(rbind, kw_results)
write.csv(kw_df, "alpha_diversity_kruskal_wallis.csv", row.names = FALSE)
message("Kruskal-Wallis results saved to: alpha_diversity_kruskal_wallis.csv")
print(kw_df)

# Beta Diversity Matrice
# Relative-abundance matrix (samples x species), already in % from MetaPhlAn
abund_relab <- t(otu_table(ps)) / 100   # convert % to proportions

# Bray-Curtis dissimilarity
bray_dist <- vegdist(abund_relab, method = "bray")

# UniFrac distances for beta diversity
if (!is.null(phy_tree(ps, errorIfNULL = FALSE))) {
  unifrac_results <- GUniFrac(abund_relab, phy_tree(ps), alpha = c(0.0, 0.5, 1.0))$unifracs
  unweighted_unifrac <- as.dist(unifrac_results[,, "d_UW"])
  weighted_unifrac   <- as.dist(unifrac_results[,, "d_1"])
  message("UniFrac distances calculated.")
} else {
  message("Skipping UniFrac (no tree available).")
  unweighted_unifrac <- NULL
  weighted_unifrac   <- NULL
}

# PCoA Ordination
run_pcoa_plot <- function(dist_obj, title_label, group_col, meta_df, filename) {
  pcoa_res  <- cmdscale(dist_obj, eig = TRUE, k = 2)
  var_exp   <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)
  pcoa_df   <- data.frame(
    PC1     = pcoa_res$points[, 1],
    PC2     = pcoa_res$points[, 2],
    SampleID = rownames(pcoa_res$points)
  )
  pcoa_df <- merge(pcoa_df, meta_df, by.x = "SampleID", by.y = sample_id_col)
  
  p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
    geom_point(size = 4, alpha = 0.85) +
    stat_ellipse(aes(group = .data[[group_col]]), level = 0.95, linetype = 2) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title_label,
      x     = paste0("PC1 (", var_exp[1], "%)"),
      y     = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right")
  
  ggsave(filename, p, width = 8, height = 6)
  invisible(list(plot = p, coords = pcoa_df))
}
 
run_pcoa_plot(bray_dist, "PCoA – Bray-Curtis",
              group_var, metadata, "pcoa_braycurtis.pdf")
 
if (!is.null(unweighted_unifrac)) {
  run_pcoa_plot(unweighted_unifrac, "PCoA – Unweighted UniFrac",
                group_var, metadata, "pcoa_unweighted_unifrac.pdf")
  run_pcoa_plot(weighted_unifrac,   "PCoA – Weighted UniFrac",
                group_var, metadata, "pcoa_weighted_unifrac.pdf")
}

# PERMANOVA (adonis2)
run_permanova <- function(dist_obj, meta_df, group_col, label) {
  formula <- as.formula(paste("dist_obj ~", group_col))
  result  <- adonis2(formula, data = meta_df, permutations = 999)
  cat("\nPERMANOVA –", label, "\n")
  print(result)
  result
}
 
perm_bray <- run_permanova(bray_dist, data.frame(metadata), group_var, "Bray-Curtis")
 
if (!is.null(unweighted_unifrac)) {
  perm_uu <- run_permanova(unweighted_unifrac, data.frame(metadata), group_var, "Unweighted UniFrac")
  perm_wu <- run_permanova(weighted_unifrac,   data.frame(metadata), group_var, "Weighted UniFrac")
}

# PERMDISP (beta-dispersion)
run_permdisp <- function(dist_obj, meta_df, group_col, label) {
  groups  <- factor(meta_df[[group_col]])
  disp    <- betadisper(dist_obj, groups)
  test    <- permutest(disp, permutations = 999)
  cat("\nPERMDISP –", label, "\n")
  print(test)
  invisible(list(disp = disp, test = test))
}
 
disp_bray <- run_permdisp(bray_dist, data.frame(metadata), group_var, "Bray-Curtis")
 
if (!is.null(unweighted_unifrac)) {
  disp_uu <- run_permdisp(unweighted_unifrac, data.frame(metadata), group_var, "Unweighted UniFrac")
  disp_wu <- run_permdisp(weighted_unifrac,   data.frame(metadata), group_var, "Weighted UniFrac")
}

# Pairwise PERMANOVA 
if (requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  pw_bray <- pairwiseAdonis::pairwise.adonis(
    bray_dist,
    factors      = data.frame(metadata)[[group_var]],
    p.adjust.m   = "BH",
    perm         = 999
  )
  print(pw_bray)
  write.csv(pw_bray, "pairwise_permanova_braycurtis.csv", row.names = FALSE)
}

write.csv(as.matrix(bray_dist), "beta_braycurtis_matrix.csv")
if (!is.null(unweighted_unifrac)) {
  write.csv(as.matrix(unweighted_unifrac), "beta_unweighted_unifrac_matrix.csv")
  write.csv(as.matrix(weighted_unifrac),   "beta_weighted_unifrac_matrix.csv")
}
