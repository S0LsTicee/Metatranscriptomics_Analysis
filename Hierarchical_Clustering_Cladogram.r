# This code utilizes Bray-Curtis and Aitchison to calculate the average and the distance between the data points 
#Install packages
my_packages <- c"tidyverse", "vegan", "ggtree", "ape",
                   "RColorBrewer", "patchwork", "ggplot2",
                   "dendextend", "factoextra", "scales"

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(vegan)
library(ggtree)
library(ape)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
library(dendextend)
library(factoextra)
library(scales)

#Load Files
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

# Auto-detect group variable
group_col <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) & n_distinct(.) > 1 & n_distinct(.) < 20)) %>%
  colnames() %>% first()

# Build Species-Abundance Matrix
colnames(metaphlan)[1] <- "taxonomy"
species_rows <- metaphlan %>%
  filter(grepl("s__", taxonomy) & !grepl("t__", taxonomy)) %>%
  mutate(species = gsub(".*s__", "s__", taxonomy))
 
sample_cols <- intersect(colnames(species_rows), metadata[[meta_id_col]])
 
abund_mat <- species_rows %>%
  select(species, all_of(sample_cols)) %>%
  column_to_rownames("species") %>%
  t() %>% as.data.frame()   # samples × species
 
# Prevalence filter: keep species in ≥ 10% of samples
prev_keep <- colSums(abund_mat > 0) / nrow(abund_mat) >= 0.10
abund_mat <- abund_mat[, prev_keep]
message(sprintf("Matrix: %d samples × %d species (prevalence ≥ 10%%)",
                nrow(abund_mat), ncol(abund_mat)))

# Distance Matrices
# Bray-Curtis on raw relative abundances
bray_dist <- vegdist(abund_mat, method = "bray")
 
# Aitchison (CLR transform + Euclidean)
clr_transform <- function(x) {
  x[x == 0] <- 0.5   # Laplace smoothing
  log(x) - mean(log(x))
}
abund_clr     <- t(apply(abund_mat, 1, clr_transform))
aitchison_dist <- dist(abund_clr, method = "euclidean")

# Hierarchical Clustering
hc_bray     <- hclust(bray_dist,     method = "ward.D2")
hc_aitchison <- hclust(aitchison_dist, method = "ward.D2")

# Convert to phylo for ggtree
phy_bray      <- as.phylo(hc_bray)
phy_aitchison <- as.phylo(hc_aitchison)

# Metadata Annotation
meta_sub <- metadata %>%
  filter(.data[[meta_id_col]] %in% sample_cols) %>%
  rename(label = all_of(meta_id_col))

# # Colour palette for group
if (!is.null(group_col) && !is.na(group_col)) {
  grps    <- unique(meta_sub[[group_col]])
  n_grps  <- length(grps)
  grp_pal <- setNames(
    colorRampPalette(brewer.pal(min(n_grps, 8), "Dark2"))(n_grps),
    grps
  )
} else {
  grp_pal <- NULL
}
 
# Detect numeric metadata columns for heatmap strips
numeric_cols <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(is.numeric)) %>%
  colnames()

# Build Annotated ggtree Plot
make_hclust_plot <- function(phy, title_str, layout = "rectangular") {
  p <- ggtree(phy, layout = layout, branch.length = "branch.length",
               color = "grey40", linewidth = 0.45) %<+% meta_sub
 
  if (!is.null(group_col) && !is.na(group_col)) {
    p <- p +
      geom_tippoint(aes(color = .data[[group_col]]), size = 3, alpha = 0.9) +
      geom_tiplab(size = 2.4, offset = 0.005) +
      scale_color_manual(values = grp_pal, name = group_col)
  } else {
    p <- p +
      geom_tippoint(size = 2.5, color = "#1565C0", alpha = 0.85) +
      geom_tiplab(size = 2.4, offset = 0.005)
  }
 
  p +
    labs(title = title_str,
         subtitle = sprintf("%d samples  |  Ward.D2 linkage", nrow(abund_mat))) +
    theme_tree2() +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9),
          legend.position = "right")
}

# Circular Cladogram
p_bray_circ <- make_hclust_plot(phy_bray,
  "Hierarchical Clustering  –  Bray-Curtis  (Circular)", layout = "circular")
p_ait_circ  <- make_hclust_plot(phy_aitchison,
  "Hierarchical Clustering  –  Aitchison  (Circular)",  layout = "circular")
 
ggsave("hclust_braycurtis_circular.pdf",  p_bray_circ,
       width = 14, height = 14, dpi = 300, limitsize = FALSE)
ggsave("hclust_aitchison_circular.pdf",   p_ait_circ,
       width = 14, height = 14, dpi = 300, limitsize = FALSE)
message("Saved: hclust_braycurtis_circular.pdf")

# Rectangular Cladogram
p_bray_rect <- make_hclust_plot(phy_bray,
  "Hierarchical Clustering  –  Bray-Curtis  (Rectangular)")
p_ait_rect  <- make_hclust_plot(phy_aitchison,
  "Hierarchical Clustering  –  Aitchison  (Rectangular)")
 
sample_h <- max(8, nrow(abund_mat) * 0.22)
ggsave("hclust_braycurtis_rectangular.pdf",  p_bray_rect,
       width = 12, height = sample_h, dpi = 300, limitsize = FALSE)
ggsave("hclust_aitchison_rectangular.pdf",   p_ait_rect,
       width = 12, height = sample_h, dpi = 300, limitsize = FALSE)
message("Saved: hclust_braycurtis_rectangular.pdf")

# Dendrogram With Colored Branches
make_dend_plot <- function(hc, dist_name) {
  dend <- as.dendrogram(hc)
  if (!is.null(group_col) && !is.na(group_col)) {
    sample_order <- hc$labels
    grp_vec      <- meta_sub[[group_col]][match(sample_order, meta_sub$label)]
    col_vec      <- grp_pal[grp_vec]
    labels_colors(dend) <- col_vec[order.dendrogram(dend)]
  }
  pdf(sprintf("hclust_%s_dendrogram_colored.pdf", tolower(gsub("-","", dist_name))),
      width = max(12, nrow(abund_mat) * 0.25), height = 8)
  par(mar = c(12, 4, 3, 1))
  plot(dend, main = sprintf("Hierarchical Clustering – %s (Ward.D2)", dist_name),
       ylab = "Height", cex = 0.7)
  if (!is.null(group_col) && !is.na(group_col)) {
    legend("topright", legend = names(grp_pal), fill = grp_pal,
           title = group_col, bty = "n", cex = 0.8)
  }
  dev.off()
  message(sprintf("Saved: hclust_%s_dendrogram_colored.pdf",
                  tolower(gsub("-","", dist_name))))
}
 
make_dend_plot(hc_bray,      "Bray-Curtis")
make_dend_plot(hc_aitchison, "Aitchison")

# Cophenetic Correlation Goodness of Fit
coph_bray <- cor(bray_dist,      cophenetic(hc_bray))
coph_ait  <- cor(aitchison_dist, cophenetic(hc_aitchison))
message(sprintf("Cophenetic correlation  –  Bray-Curtis: %.3f  |  Aitchison: %.3f",
                coph_bray, coph_ait))

# Silhouette Scores (if group col is detected)
if (!is.null(group_col) && !is.na(group_col)) {
  library(cluster)
  grp_int      <- as.integer(factor(meta_sub[[group_col]][match(hc_bray$labels, meta_sub$label)]))
  sil_bray     <- silhouette(grp_int, bray_dist)
  sil_ait      <- silhouette(grp_int, aitchison_dist)
  mean_sil_b   <- mean(sil_bray[, "sil_width"])
  mean_sil_a   <- mean(sil_ait[,  "sil_width"])
  message(sprintf("Mean silhouette  –  Bray-Curtis: %.3f  |  Aitchison: %.3f",
                  mean_sil_b, mean_sil_a))
 
  sil_df <- rbind(
    data.frame(sample = rownames(sil_bray), sil_width = sil_bray[,"sil_width"],
               cluster = sil_bray[,"cluster"], distance = "Bray-Curtis"),
    data.frame(sample = rownames(sil_ait),  sil_width = sil_ait[,"sil_width"],
               cluster = sil_ait[,"cluster"],  distance = "Aitchison")
  )
 
  p_sil <- ggplot(sil_df, aes(x = reorder(sample, sil_width), y = sil_width,
                               fill = factor(cluster))) +
    geom_col() +
    facet_wrap(~ distance, scales = "free_x") +
    scale_fill_brewer(palette = "Dark2", name = "Cluster") +
    coord_flip() +
    labs(title = "Silhouette Plot – Cluster Quality",
         x = NULL, y = "Silhouette Width") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
 
  ggsave("hclust_silhouette.pdf", p_sil,
         width = 12, height = max(6, nrow(abund_mat) * 0.15), dpi = 300)
  message("Saved: hclust_silhouette.pdf")
}

# Export distance matrices & cluster assignments 
# Optimal k by silhouette (k = 2 to min(10, n-1))
fviz_nbclust_safe <- function(mat, max_k = min(10, nrow(mat)-1)) {
  tryCatch({
    wss <- sapply(2:max_k, function(k) {
      km <- kmeans(mat, centers = k, nstart = 25)
      km$tot.withinss
    })
    data.frame(k = 2:max_k, wss = wss)
  }, error = function(e) NULL)
}
 
wss_df <- fviz_nbclust_safe(abund_clr)
if (!is.null(wss_df)) {
  p_elbow <- ggplot(wss_df, aes(x = k, y = wss)) +
    geom_line(color = "#1565C0") + geom_point(size = 3, color = "#1565C0") +
    labs(title = "Elbow Plot – Optimal Number of Clusters",
         x = "Number of clusters (k)", y = "Total within-cluster SS") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  ggsave("hclust_elbow_plot.pdf", p_elbow, width = 7, height = 5, dpi = 300)
  message("Saved: hclust_elbow_plot.pdf")
}
 
# Cut tree at k = 2 (default) and export
k_cut    <- 2
clusters_bray <- cutree(hc_bray, k = k_cut)
clusters_ait  <- cutree(hc_aitchison, k = k_cut)
 
cluster_df <- data.frame(
  sample_id       = names(clusters_bray),
  cluster_bray    = clusters_bray,
  cluster_aitchison = clusters_ait
) %>%
  left_join(meta_sub %>% rename(sample_id = label), by = "sample_id")
 
write.csv(cluster_df, "hclust_cluster_assignments.csv", row.names = FALSE)
write.csv(as.matrix(bray_dist),      "distance_braycurtis.csv")
write.csv(as.matrix(aitchison_dist), "distance_aitchison.csv")
 
message("Saved: All Matrices")
