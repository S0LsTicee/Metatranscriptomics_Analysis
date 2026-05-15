#Install Packages
my_packages <- c("tidyverse", "vegan", "dendextend", "ggdendro",
                   "ggplot2", "patchwork", "RColorBrewer",
                   "cluster", "factoextra", "scales", "circlize",
                   "ComplexHeatmap")

install.packages(my_packages)

# Check whether the required packages are installed and install them if not yet installed
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library(tidyverse)
library(vegan)
library(dendextend)
library(ggdendro)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(cluster)
library(factoextra)
library(scales)
library(circlize)
library(ComplexHeatmap)

# Load Files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder of your merged metaphlan file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE, 
)

metadata <- read.delim(
  "metadata_merged.txt",  #This is a placeholder of your merged metadata file
  header = TRUE,
  sep = "\t", 
  check.names = FALSE
)

meta_id_col   <- colnames(metadata)[1]

#Build Species Abundance Matrix
colnames(metaphlan)[1] <- "taxonomy"
 
abund_mat <- metaphlan %>%
  filter(grepl("s__", taxonomy) & !grepl("t__", taxonomy)) %>%
  mutate(species = gsub(".*s__", "s__", taxonomy)) %>%
  select(species, any_of(metadata[[meta_id_col]])) %>%
  column_to_rownames("species") %>%
  t() %>% as.data.frame()
 
# Prevalence filter ≥ 10 %
prev_keep <- colSums(abund_mat > 0) / nrow(abund_mat) >= 0.10
abund_mat <- abund_mat[, prev_keep]
message(sprintf("Matrix: %d samples × %d species", nrow(abund_mat), ncol(abund_mat)))
 
# Align metadata
meta_aligned <- metadata[match(rownames(abund_mat), metadata[[meta_id_col]]), ]

# Auto-detect Annotation Columns
cat_cols <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(~ !is.numeric(.) & n_distinct(.) > 1 & n_distinct(.) < 20)) %>%
  colnames()
 
num_cols <- metadata %>%
  select(-all_of(meta_id_col)) %>%
  select(where(is.numeric)) %>%
  colnames()
 
group_col <- cat_cols[1]   # primary grouping variable

# Distance Matrices
# Bray-Curtis
bray_dist <- vegdist(abund_mat, method = "bray")
 
# Jaccard (binary)
jacc_dist <- vegdist(abund_mat, method = "jaccard", binary = TRUE)
 
# Aitchison (CLR + Euclidean)
clr_mat <- t(apply(abund_mat, 1, function(x) {
  x[x == 0] <- 0.5
  log(x) - mean(log(x))
}))
ait_dist <- dist(clr_mat, method = "euclidean")

# Color Palettes
build_pal <- function(col_name, df = meta_aligned) {
  if (is.null(col_name) || is.na(col_name)) return(NULL)
  vals <- unique(df[[col_name]])
  n    <- length(vals)
  pal  <- setNames(
    colorRampPalette(brewer.pal(min(n, 8), "Dark2"))(n),
    vals
  )
  pal
}
 
grp_pal <- build_pal(group_col)

# Color Labels by Group
label_colors <- if (!is.null(grp_pal)) {
  grp_pal[meta_aligned[[group_col]][match(rownames(abund_mat),
                                          meta_aligned[[meta_id_col]])]]
} else {
  rep("grey30", nrow(abund_mat))
}
names(label_colors) <- rownames(abund_mat)

# Core Dendrogram Function
save_dendrogram <- function(dist_obj, method = "ward.D2",
                             dist_name = "BrayCurtis",
                             link_name = "Ward.D2") {
 
  hc   <- hclust(dist_obj, method = method)
  dend <- as.dendrogram(hc)
 
  # Colour branches by group (requires dendextend)
  if (!is.null(grp_pal) && !is.na(group_col)) {
    grp_vec <- meta_aligned[[group_col]][match(hc$labels, meta_aligned[[meta_id_col]])]
    dend    <- color_branches(dend, k = length(unique(grp_vec)),
                               col = grp_pal[unique(grp_vec)])
    labels_colors(dend) <- label_colors[labels(dend)]
  }
 
  # Cophenetic correlation
  coph <- round(cor(dist_obj, cophenetic(hc)), 3)
 
  fname <- sprintf("dendrogram_%s_%s.pdf", dist_name, link_name)
  pdf(fname, width = max(12, nrow(abund_mat) * 0.3), height = 9)
  par(mar = c(12, 5, 4, 2))
  plot(dend,
       main  = sprintf("Sample Clustering  –  %s  |  %s linkage", dist_name, link_name),
       sub   = sprintf("Cophenetic r = %s  |  n = %d samples", coph, nrow(abund_mat)),
       ylab  = "Distance",
       cex   = 0.75,
       xlab  = "")
  if (!is.null(grp_pal) && !is.na(group_col)) {
    legend("topright", legend = names(grp_pal), fill = grp_pal,
           title = group_col, bty = "n", cex = 0.75)
  }
  dev.off()
  message(sprintf("Saved: %s  (cophenetic r = %s)", fname, coph))
 
  list(hc = hc, dend = dend, coph = coph)
}

# Generate Dendrograms All Distance-Linkage Combos
distances <- list(
  list(dist = bray_dist, name = "BrayCurtis"),
  list(dist = jacc_dist, name = "Jaccard"),
  list(dist = ait_dist,  name = "Aitchison")
)
 
linkages <- c("ward.D2", "average", "complete")
 
results <- list()
for (d in distances) {
  for (lnk in linkages) {
    key          <- paste0(d$name, "_", gsub("\\.", "", lnk))
    lnk_label    <- tools::toTitleCase(gsub("\\.", " ", lnk))
    results[[key]] <- save_dendrogram(d$dist, method = lnk,
                                       dist_name = d$name,
                                       link_name = lnk_label)
  }
}

# ggplot2 Dendrogram
# Best dendrogram = highest cophenetic r among Ward.D2 results
coph_vals <- sapply(c("BrayCurtis_wardD2", "Jaccard_wardD2", "Aitchison_wardD2"),
                    function(k) results[[k]]$coph)
best_key  <- names(which.max(coph_vals))
best_hc   <- results[[best_key]]$hc
best_name <- gsub("_wardD2", "", best_key)
message(sprintf("Best dendrogram (highest cophenetic r): %s", best_key))
 
dend_data  <- dendro_data(best_hc, type = "rectangle")
label_df   <- dend_data$labels %>%
  left_join(meta_aligned %>% rename(label = all_of(meta_id_col)), by = "label")
 
p_ggdend <- ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "grey40", linewidth = 0.5) +
  {if (!is.null(group_col) && !is.na(group_col))
    geom_text(data = label_df,
              aes(x = x, y = y - max(dend_data$segments$y) * 0.03,
                  label = label, color = .data[[group_col]]),
              angle = 90, hjust = 1, size = 2.8)
  else
    geom_text(data = label_df,
              aes(x = x, y = y - max(dend_data$segments$y) * 0.03, label = label),
              angle = 90, hjust = 1, size = 2.8, color = "grey30")
  } +
  {if (!is.null(grp_pal))
    scale_color_manual(values = grp_pal, name = group_col)
  } +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.05))) +
  labs(
    title    = sprintf("Sample Clustering Dendrogram  –  %s  |  Ward.D2", best_name),
    subtitle = sprintf("Cophenetic r = %s  |  %d samples  |  %d species",
                       coph_vals[best_key], nrow(abund_mat), ncol(abund_mat)),
    x        = NULL,
    y        = "Distance"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9),
    legend.position  = "right"
  )
 
ggsave("dendrogram_best_ggplot.pdf", p_ggdend,
       width = max(10, nrow(abund_mat) * 0.3), height = 7, dpi = 300)
message("Saved: dendrogram_best_ggplot.pdf")

#ComplexHeatmap with Dendrogram for Top 40 Most Variable Species
top_sp <- names(sort(apply(abund_mat, 2, var), decreasing = TRUE)[1:min(40, ncol(abund_mat))])
heat_mat <- t(as.matrix(abund_mat[, top_sp]))
heat_mat_scale <- t(scale(t(heat_mat)))   # z-score per species
 
# Column (sample) annotation
if (!is.null(group_col) && !is.na(group_col)) {
  ann_df  <- meta_aligned %>%
    select(all_of(meta_id_col), all_of(group_col)) %>%
    column_to_rownames(meta_id_col)
  col_ann <- HeatmapAnnotation(
    df  = ann_df,
    col = setNames(list(grp_pal), group_col),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 9)
  )
} else {
  col_ann <- NULL
}
 
heat_col <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#D6604D")
)
 
ht <- Heatmap(
  heat_mat_scale,
  name                  = "Z-score",
  col                   = heat_col,
  top_annotation        = col_ann,
  cluster_columns       = best_hc,
  clustering_distance_rows = "euclidean",
  clustering_method_rows   = "ward.D2",
  show_row_names        = TRUE,
  show_column_names     = TRUE,
  row_names_gp          = gpar(fontsize = 7),
  column_names_gp       = gpar(fontsize = 8),
  column_names_rot      = 45,
  row_names_max_width   = unit(6, "cm"),
  heatmap_legend_param  = list(title = "Z-score", title_gp = gpar(fontsize = 9)),
  column_title          = sprintf("Top %d variable species  |  Samples clustered by %s",
                                  length(top_sp), best_name),
  column_title_gp       = gpar(fontsize = 10, fontface = "bold")
)
 
pdf("dendrogram_heatmap.pdf",
    width  = max(12, nrow(abund_mat) * 0.35),
    height = max(10, length(top_sp) * 0.25))
draw(ht)
dev.off()
message("Saved: dendrogram_heatmap.pdf")

# Sihouette & Elbow Plots 
# Silhouette for k = 2 to min(10, n-1)
max_k   <- min(10, nrow(abund_mat) - 1)
best_dist_obj <- switch(best_name,
  BrayCurtis = bray_dist, Jaccard = jacc_dist, Aitchison = ait_dist,
  bray_dist)
 
sil_scores <- sapply(2:max_k, function(k) {
  cl  <- cutree(best_hc, k = k)
  sil <- silhouette(cl, best_dist_obj)
  mean(sil[, "sil_width"])
})
 
p_sil <- ggplot(data.frame(k = 2:max_k, sil = sil_scores),
                aes(x = k, y = sil)) +
  geom_line(color = "#E53935", linewidth = 0.9) +
  geom_point(size = 3, color = "#E53935") +
  geom_vline(xintercept = which.max(sil_scores) + 1,
             linetype = "dashed", color = "grey40") +
  labs(title = "Optimal Cluster Number  –  Silhouette Width",
       x = "Number of clusters (k)", y = "Mean silhouette width") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
 
ggsave("dendrogram_silhouette.pdf", p_sil, width = 7, height = 5, dpi = 300)
message(sprintf("Saved: dendrogram_silhouette.pdf  (optimal k = %d)",
                which.max(sil_scores) + 1))

# Cophenetic Correlation Summary
                    coph_summary <- data.frame(
  distance = rep(c("BrayCurtis","Jaccard","Aitchison"), each = 3),
  linkage  = rep(c("Ward.D2","Average","Complete"), 3),
  cophenetic_r = sapply(names(results), function(k) results[[k]]$coph)
)
print(coph_summary %>% arrange(desc(cophenetic_r)))
write.csv(coph_summary, "dendrogram_cophenetic_summary.csv", row.names = FALSE)
message("Saved: dendrogram_cophenetic_summary.csv")

# Cluster Assignment at Optimal k
opt_k      <- which.max(sil_scores) + 1
cluster_df <- data.frame(
  sample_id = rownames(abund_mat),
  cluster   = cutree(best_hc, k = opt_k)
) %>% left_join(meta_aligned %>% rename(sample_id = all_of(meta_id_col)),
                by = "sample_id")
 
write.csv(cluster_df, "dendrogram_cluster_assignments.csv", row.names = FALSE)
message(sprintf("Saved: dendrogram_cluster_assignments.csv  (k = %d)", opt_k))
