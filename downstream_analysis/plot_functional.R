library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(pheatmap)
library(vegan)

meta <- read.csv("metadata.csv", stringsAsFactors = FALSE) %>%
  mutate(group = factor(group, levels = c("Healthy", "Crohn")))
meta <- meta[meta$sample != "SRR5983484", ]

# tabella HUMAnN merged
hum <- read.delim("humann/pathabundance_merged_cpm.tsv",
                  check.names = FALSE,
                  stringsAsFactors = FALSE)

names(hum)[1] <- "feature"
names(hum)[-1] <- sub("\\.clean\\.concat_Abundance-CPM$", "", names(hum)[-1])


hum_long <- hum %>%
  pivot_longer(-feature, names_to = "sample", values_to = "abundance") %>%
  mutate(
    is_stratified = str_detect(feature, fixed("|")),
    pathway = sub("\\|.*$", "", feature),
    taxon = ifelse(is_stratified, sub("^.*\\|", "", feature), NA)
  ) %>%
  left_join(meta, by = "sample")


special_rows <- c("UNMAPPED", "UNINTEGRATED")

# unstratified pathway
path_unstrat <- hum_long %>%
  filter(!is_stratified, !pathway %in% special_rows)

# stratified pathway
path_strat <- hum_long %>%
  filter(is_stratified, !pathway %in% special_rows)


mat_path_df <- path_unstrat %>%
  select(pathway, sample, abundance) %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)

mat_path <- as.data.frame(mat_path_df)
rownames(mat_path) <- mat_path$pathway
mat_path$pathway <- NULL
mat_path <- as.matrix(mat_path)

grp <- meta$group[match(colnames(mat_path), meta$sample)]

#PCA
X <- t(mat_path)
X_log <- log10(X + 1)
X_log <- X_log[, apply(X_log, 2, var) > 0, drop = FALSE]

pca <- prcomp(X_log, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_df$group <- meta$group[match(pca_df$sample, meta$sample)]

pca_df

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample), vjust = -0.7, size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(
    title = "PCA of HUMAnN pathway",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
  )


#Preprocessing Heatmap
pc <- 1               # pseudocount per log10
prev_cutoff <- 1     # una pathway è "presente" se >= 10 CPM
topN <- 25

mean_crohn <- rowMeans(mat_path[, grp == "Crohn", drop = FALSE], na.rm = TRUE)
mean_ctrl  <- rowMeans(mat_path[, grp == "Healthy", drop = FALSE], na.rm = TRUE)

log2FC <- log2((mean_crohn + 1) / (mean_ctrl + 1))
direction <- ifelse(log2FC > 0, "Enriched in Crohn", "Depleted in Crohn")

pval <- apply(mat_path, 1, function(x) {
  wilcox.test(
    x[grp == "Crohn"],
    x[grp == "Healthy"],
    exact = FALSE
  )$p.value
})

padj <- p.adjust(pval, method = "BH")

prev_crohn <- rowSums(mat_path[, grp == "Crohn", drop = FALSE] >= prev_cutoff)
prev_ctrl  <- rowSums(mat_path[, grp == "Healthy", drop = FALSE] >= prev_cutoff)

abs_log2FC <- abs(log2FC)

res_path <- tibble(
  pathway = rownames(mat_path),
  log2FC = log2FC,
  direction = direction,
  pval = pval,
  padj = padj,
  prev_crohn = prev_crohn,
  prev_ctrl = prev_ctrl,
  abs_log2FC = abs_log2FC
) %>%
  arrange(padj, desc(abs_log2FC))

print(head(res_path, 20))

# Global top 50

topN_global <- 50

top_pathways <- res_path %>%
  filter(prev_crohn >= 2 & prev_ctrl >= 2) %>%
  slice_head(n = topN_global) %>%
  pull(pathway)

mat_plot <- log10(mat_path[top_pathways, , drop = FALSE] + pc)

ann_col <- data.frame(Group = grp)
rownames(ann_col) <- colnames(mat_plot)

pheatmap(
  mat_plot,
  border_color = "black",
  annotation_col = ann_col,
  main = paste0(
      "                                                                                   ",
      "Top ", topN_global," pathways Crohn vs Healthy (prevalence >= 2 in both groups)"
  ),
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  cellwidth = 20,
  cellheight = 10
)


# Top 25 enriched in Crohn
enriched_pathways <- res_path %>%
  filter(log2FC > 0, prev_crohn >= 1) %>%
  arrange(desc(log2FC), desc(prev_crohn)) %>%
  slice_head(n = topN) %>%
  pull(pathway)

# Top 25 enriched in Healthy

depleted_pathways <- res_path %>%
  filter(log2FC < 0, prev_ctrl >= 1) %>%
  arrange(desc(abs_log2FC), desc(prev_ctrl)) %>%
  slice_head(n = topN) %>%
  pull(pathway)


#Heatmap enriched/depleted

mat_enriched <- log10(mat_path[enriched_pathways, , drop = FALSE] + pc)
mat_depleted <- log10(mat_path[depleted_pathways, , drop = FALSE] + pc)

ann_row <- data.frame(Group = grp)
rownames(ann_row) <- colnames(mat_enriched)

pheatmap(
  (mat_enriched),
  annotation_col = ann_row,
  main = paste0(
    "                                                                                   ",
    "Top ", topN," pathways enriched in Crohn (prevalence >= 1 in Crohn)"
  ),
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 315,
  border_color = "black",
  cellwidth = 20,
  cellheight = 10
)

pheatmap(
  (mat_depleted),
  annotation_col = ann_row,
  main = paste0(
    "                                                                                   ",
    "Top ", topN, " pathways enriched in Healthy (prevalence >= 1 in Healthy)"
  ),
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 315,
  border_color = "black",
  cellwidth = 20,
  cellheight = 10
)





# Contributor for pathway PWY4140
target_pathway <- top_pathways[16]
path_total <- hum_long %>%
  filter(pathway == target_pathway, !is_stratified) %>%
  select(sample, total_pathway = abundance)

# Txa contribution
contrib <- hum_long %>%
  filter(pathway == target_pathway, is_stratified) %>%
  left_join(path_total, by = "sample") %>%
  mutate(
    pct_of_pathway = 100 * abundance / total_pathway
  )

top_taxa <- contrib %>%
  filter(taxon != "unclassified") %>%
  group_by(taxon) %>%
  summarise(mean_pct = mean(pct_of_pathway, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_pct)) %>%
  #slice_head(n = 8) %>%
  pull(taxon)

contrib_plot <- contrib %>%
  mutate(
    taxon_plot = case_when(
      taxon == "unclassified" ~ "unclassified",
      taxon %in% top_taxa ~ taxon,
      TRUE ~ "Other"
    )
  ) %>%
  group_by(sample, group, taxon_plot) %>%
  summarise(pct_of_pathway = sum(pct_of_pathway), .groups = "drop")

# stacked bar per sample
ggplot(contrib_plot, aes(x = sample, y = pct_of_pathway, fill = taxon_plot)) +
  geom_col() +
  facet_wrap(~group, scales = "free_x") +
  coord_flip() +
  theme_bw() +
  labs(title = paste0("Taxonomic contributors (", target_pathway,")"),
       y = "% of pathway",
       x = "Sample")

# stacked bar per group
contrib_group <- contrib_plot %>%
  group_by(group, taxon_plot) %>%
  summarise(mean_pct = mean(pct_of_pathway), .groups = "drop")

ggplot(contrib_group, aes(x = group, y = mean_pct, fill = taxon_plot)) +
  geom_col() +
  theme_bw() +
  labs(title = paste0("Mean contributors by group (", target_pathway,")"),
       y = "Mean % of pathway",
       x = "")


