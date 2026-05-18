#Install packages
packages <- c("ggplot2", "dplyr", "tidyr", "tibble",
          "vegan", "patchwork", "RColorBrewer", "ggrepel")

install.packages(packages)

# Check for missing files and install them
not_installed <- my_packages[!(my_packages %in% installed.packages()[, "Package"])]

if(length(not_installed)) {
  install.packages(not_installed)
}

library("ggplot2")
library("dplr")
library("tidyr")
library("tibble")
library("vegan")
library("patchwork")
library("RColorBrewer")
library("ggrepel")

# Load merged files
metaphlan <- read.delim(
  "metaphlan_merged.txt", #This is a placeholder for the name of your merged metaphlan file
  header           = TRUE,
  sep              = "\t",
  comment.char     = "#",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.delim(
  "metadata_merged.txt", #This is a placeholder for the name of your merged metaphlan file
  header           = TRUE,
  sep              = "\t",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

#Parse metadata to avoid column count error
metadata_p <- function(path) {
raw_lines <- readLines(path, warn = FALSE)
  
# Remove BOM if present
raw_lines[1] <- gsub("^\xef\xbb\xbf", "", raw_lines[1])
split_lines <- strsplit(raw_lines, ",")

data_rows <- split_lines[-1]   # drop header line

df <- data.frame(
    sample_id   = sapply(data_rows, function(r) r[1]),
    Sample.Name = sapply(data_rows, function(r) r[length(r) - 1]),
    SRA.Study   = sapply(data_rows, function(r) r[length(r)]),
    stringsAsFactors = FALSE
  )
df$group <- sub("_\\d+$", "", df$Sample.Name)   # Drops "_" in column name
  df$group <- gsub("_", " ", df$group)             # Drops "_" in rows
  df
}

meta <- read_metadata_safe("master_metadata.csv")
stopifnot(all(grepl("^SRR", meta$sample_id))) 
message("Metadata loaded: ", nrow(meta), " samples | groups: ",
        paste(sort(unique(meta$group)), collapse = ", "))

# Strip "_profile.txt" suffix from sample column names
colnames(metaphlan) <- gsub("_profile\\.txt$", "", colnames(metaphlan))

# Keep species level only (contains s__, excludes t__ strain level)
sp <- metaphlan[
  grepl("s__", mphcommon <- intersect(colnames(sp_bin), meta$sample_id)
if (length(common) == 0)
  stop("No matching sample IDs between MetaPhlAn and metadata.\n",
       "  MetaPhlAn cols: ", paste(head(colnames(sp_bin), 5), collapse=", "), "\n",
       "  Metadata IDs:   ", paste(head(meta$sample_id,   5), collapse=", "))
 
sp_bin <- sp_bin[, common, drop = FALSE]
meta   <- meta[match(common, meta$sample_id), ]
rownames(meta) <- NULL_raw$clade_name) & !grepl("t__", mph_raw$clade_name),
]
rownames(sp) <- sp$clade_name
sp$clade_name <- NULL

# Convert 0/2 encoding to binary 0/1
sp_bin <- as.data.frame(lapply(sp, function(x) as.integer(x > 0)))
rownames(sp_bin) <- rownames(sp)

#Align Samples
common <- intersect(colnames(sp_bin), meta$sample_id)
if (length(common) == 0)
  stop("No matching sample IDs between MetaPhlAn and metadata.\n",
       "  MetaPhlAn cols: ", paste(head(colnames(sp_bin), 5), collapse=", "), "\n",
       "  Metadata IDs:   ", paste(head(meta$sample_id,   5), collapse=", "))
 
sp_bin <- sp_bin[, common, drop = FALSE]
meta   <- meta[match(common, meta$sample_id), ]
rownames(meta) <- NULL

#Group vectors & Color palette
GROUP_COL <- "group"   # change if your grouping column differs
groups    <- sort(unique(meta[[GROUP_COL]]))
pal       <- setNames(
  colorRampPalette(brewer.pal(min(max(length(groups),3),8),"Set1"))(length(groups)),
  groups
)
# Override with publication-friendly colours for the two known groups
if (setequal(groups, c("Group1","Group2")))
  pal <- c("Group1" = "#2196F3", "Group2" = "#FF5722")
 
g1_samp <- meta$sample_id[meta[[GROUP_COL]] == "Group1"]
g2_samp  <- meta$sample_id[meta[[GROUP_COL]] == "Group2"]
 
short_sp  <- function(x) sub(".*s__", "", x)   # shorten clade_name to species

#Panel A: Species Richness
rich_df <- data.frame(
  sample_id = colnames(sp_bin),
  Richness  = colSums(sp_bin),
  stringsAsFactors = FALSE
)
rich_df <- merge(rich_df, meta, by = "sample_id")
 
p_rich <- ggplot(rich_df,
    aes(x = .data[[GROUP_COL]], y = Richness,
        fill = .data[[GROUP_COL]], colour = .data[[GROUP_COL]])) +
  geom_boxplot(alpha = 0.5, width = 0.45, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = Sample.Name),
                  size = 2.8, colour = "grey30", seed = 1) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal) +
  labs(title    = "A  Species Richness",
       subtitle = "Detected species per sample",
       x = NULL, y = "Observed species") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey50"))

#Panel B: Prevalence 
prev_df <- data.frame(
  prevalence = rowSums(sp_bin),
  stringsAsFactors = FALSE
)
 
p_prev <- ggplot(prev_df, aes(x = factor(prevalence))) +
  geom_bar(fill = "#607D8B", colour = "white", width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.4, size = 3.2, colour = "grey30") +
  labs(title    = "B  Prevalence Distribution",
       subtitle = "# samples each species is detected in (max 12)",
       x = "Number of samples", y = "Number of species") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey50"))

#Panel C: Jaccard PCoA
jacc     <- vegdist(t(sp_bin), method = "jaccard", binary = TRUE)
pcoa_res <- cmdscale(jacc, k = 2, eig = TRUE)
eigs     <- pcoa_res$eig
pct_var  <- round(eigs / sum(eigs[eigs > 0]) * 100, 1)
 
pcoa_df <- data.frame(
  sample_id = colnames(sp_bin),
  PC1       = pcoa_res$points[, 1],
  PC2       = pcoa_res$points[, 2],
  stringsAsFactors = FALSE
)
pcoa_df <- merge(pcoa_df, meta, by = "sample_id")
 
p_pcoa <- ggplot(pcoa_df,
    aes(x = PC1, y = PC2,
        colour = .data[[GROUP_COL]], fill = .data[[GROUP_COL]])) +
  stat_ellipse(geom = "polygon", alpha = 0.10, type = "t",
               level = 0.90, linewidth = 0.5) +
  geom_point(size = 4, alpha = 0.95) +
  geom_text_repel(aes(label = Sample.Name),
                  size = 2.8, colour = "grey25",
                  max.overlaps = 20, seed = 42) +
  scale_colour_manual(values = pal, name = "Group") +
  scale_fill_manual(values   = pal, name = "Group") +
  labs(title    = "C  Beta Diversity – Jaccard PCoA",
       subtitle = "Binary presence/absence; 90% confidence ellipses",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey50"))

#Panel D: Shared v.s. Group-specific Species Comparison
g1_sp   <- rownames(sp_bin)[rowSums(sp_bin[, g1_samp]) > 0]
g2_sp    <- rownames(sp_bin)[rowSums(sp_bin[, g2_samp])  > 0]
shared   <- intersect(g1_sp, g2_sp)
g1_only <- setdiff(g1_sp, g2_sp)
g2_only  <- setdiff(g2_sp,  g1_sp)
 
venn_df <- data.frame(
  Category = factor(
    c("Group1 only", "Shared", "Group2 only"),
    levels = c("Group1 only", "Shared", "Group2 only")
  ),
  Count    = c(length(g1_only), length(shared), length(g2_only)),
  fill_grp = c("Group1", "Shared", "Group2")
)
pal_v <- c("Group1" = "#2196F3", "Shared" = "#9E9E9E",
           "Group2" = "#FF5722")
 
p_venn <- ggplot(venn_df, aes(x = Category, y = Count, fill = fill_grp)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(aes(label = Count), vjust = -0.4, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = pal_v) +
  labs(title    = "D  Shared vs Group-Specific Species",
       subtitle = paste0("Total: ", length(union(g1_sp, g2_sp)),
                         " unique species detected"),
       x = NULL, y = "Number of species") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey50"))

#Panel E: Top Differentially Prevalent Species
diff_df <- data.frame(
  species  = rownames(sp_bin),
  g1_prev = rowSums(sp_bin[, g1_samp]) / length(g1_samp),
  g2_prev  = rowSums(sp_bin[, g2_samp])  / length(g2_samp),
  stringsAsFactors = FALSE
)
diff_df$diff      <- diff_df$g1_prev - diff_df$g2_prev
diff_df$abs_diff  <- abs(diff_df$diff)
diff_df$direction <- ifelse(diff_df$diff > 0, "Group1", "Group2")
diff_df$label     <- short_sp(diff_df$species)
diff_df <- diff_df[order(-diff_df$abs_diff), ][1:20, ]
 
p_diff <- ggplot(diff_df,
    aes(x = diff, y = reorder(label, diff), fill = direction)) +
  geom_col(width = 0.7, colour = "white") +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
  scale_fill_manual(values = pal, name = "Higher in") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title    = "E  Top 20 Differentially Prevalent Species",
       subtitle = "Difference in detection rate between groups",
       x = "Prevalence difference (Group1 – Group2)",
       y = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y   = element_text(size = 8.5),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey50"))

#Panel F: Phylum-level Detection
phylum_df <- metaphlan
colnames(phylum_df) <- gsub("_profile\\.txt$", "", colnames(phylum_df))
phylum_df <- phylum_df[
  grepl("p__",  phylum_df$clade_name) &
  !grepl("c__", phylum_df$clade_name) &
  !grepl("UNCLASSIFIED", phylum_df$clade_name),
]
phylum_df$phylum <- sub(".*p__", "", phylum_df$clade_name)
phylum_df <- phylum_df[, c("phylum", common)]
 
phylum_long <- tidyr::pivot_longer(phylum_df,
                                   cols      = -phylum,
                                   names_to  = "sample_id",
                                   values_to = "present")
phylum_long$present <- as.integer(phylum_long$present > 0)
phylum_long <- merge(phylum_long, meta, by = "sample_id")
 
phylum_sum <- aggregate(present ~ group + phylum,
                        data = phylum_long, FUN = sum)
 
p_phylum <- ggplot(phylum_sum,
    aes(x = reorder(phylum, present),
        y = present, fill = group)) +
  geom_col(position = "dodge", width = 0.65, colour = "white") +
  scale_fill_manual(values = pal, name = "Group") +
  coord_flip() +
  labs(title    = "F  Phylum-Level Detection",
       subtitle = "Total sample-level presences per phylum",
       x = NULL, y = "Sum of sample detections") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey50"))

#Combine with patchwork
layout <- "
AABB
CCCC
DDEE
FFFF
"
 
combined <- (p_rich + p_prev + p_pcoa + p_venn + p_diff + p_phylum) +
  plot_layout(design = layout) +
  plot_annotation(
    title    = "Mouse Gut Microbiome – Taconic vs Charles River",
    subtitle = paste0(nrow(sp_bin), " species  |  ",
                      ncol(sp_bin), " samples  |  SRP075802"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey45")
    )
  )

#Save output files
ggsave("microbiome_multipanel.pdf",
       plot = combined, width = 17, height = 22, units = "in", dpi = 300)
ggsave("microbiome_multipanel.png",
       plot = combined, width = 17, height = 22, units = "in", dpi = 150)
 
message("\n Done! Saved: microbiome_multipanel.pdf / .png")
