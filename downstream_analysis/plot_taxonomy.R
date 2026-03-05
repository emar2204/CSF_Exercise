library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(vegan)
library(pheatmap)
library(tibble)
library(VennDiagram)
library(ggrepel)
#
# READ INPUT DATA
#
infile <- "metaphlan/merged_abundance_table.tsv"
mpa <- read.table(infile,  header = TRUE,sep = "\t",quote = "")
colnames(mpa) <- sub("\\.profile$", "", colnames(mpa))
meta <- read_csv("metadata.csv", show_col_types = FALSE)
#REMOVE LOW QUALITY SAMPLE
#mpa <- mpa %>% select(-`SRR5983484`)
meta <- meta[meta$sample != "SRR5983484", ]

tax_col <- names(mpa)[1]

#Extract Species only
mpa_s <- mpa %>%
  filter(str_detect(.data[[tax_col]], "\\|s__[^|]+$")) %>%
  mutate(
    species = str_extract(.data[[tax_col]], "s__[^|]+$") %>%
      str_remove("^s__") %>%
      str_replace_all("_", " ")
  )

mat <- mpa_s %>%
  select(-all_of(tax_col)) %>%
  mutate(across(-species, ~ suppressWarnings(as.numeric(.)))) %>%
  group_by(species) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("species") %>%
  as.matrix()

#Align matrix with metadata
group_vec <- meta$group
names(group_vec) <- meta$sample

keep_samp <- intersect(colnames(mat), names(group_vec))
mat2 <- mat[, keep_samp, drop = FALSE]
grp <- factor(group_vec[keep_samp], levels = c("Healthy", "Crohn"))
stopifnot(all(!is.na(grp)))


# Select those rows with at least one sample per group
present <- mat > 0
keep_taxa <- apply(present, 1, function(x) max(tapply(x, grp, sum)) >= 1)
mat_f <- mat[keep_taxa, , drop = FALSE]


#Venn Diagram
prev_crohn   <- rowSums(mat_f[, grp=="Crohn",  drop=FALSE] > 0)
prev_healthy <- rowSums(mat_f[, grp=="Healthy",drop=FALSE] > 0)

A <- rownames(mat_f)[prev_crohn > 0]
B <- rownames(mat_f)[prev_healthy > 0]

grid.newpage()
grid.draw(venn.diagram(
  x = list(Crohn = A, Healthy = B),
  filename = NULL,     # draw to current device
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.5,
  main.cex = 1.5,
  main = "Species overlap between group"
))





# Alpha diversity

alpha <- data.frame(
  sample   = colnames(mat),
  shannon  = diversity(t(mat), index = "shannon"),
  simpson  = diversity(t(mat), index = "simpson"),
  richness = specnumber(t(mat))
) %>% left_join(meta, by="sample")

p_alpha1 <- ggplot(alpha, aes(x=group, y=shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0) +
  theme_bw() + labs(title=paste0("Alpha diversity (Shannon) - Species"))
plot(p_alpha1)

p_alpha2 <- ggplot(alpha, aes(x=group, y=richness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0) +
  theme_bw() + labs(title=paste0("Alpha diversity (Richness) - Species"))
plot(p_alpha2)

alpha_tests <- list(
  shannon_wilcox  = wilcox.test(shannon ~ group, data = alpha),
  richness_wilcox = wilcox.test(richness ~ group, data = alpha)
)



# Beta diversity
bray <- vegdist(t(mat), method = "bray")

# PERMANOVA
perm <- adonis2(bray ~ group, data = meta)
perm

#PCoA
pcoa <- cmdscale(bray, k = 2, eig = TRUE)

var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100

df_pcoa <- tibble(
  sample = rownames(pcoa$points),
  PC1 = pcoa$points[, 1],
  PC2 = pcoa$points[, 2]
) %>%
  left_join(meta, by = "sample")

# Plot
ggplot(df_pcoa, aes(PC1, PC2, color = group, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  labs(
    title = "PCoA (Bray-Curtis)",
    x = sprintf("PC1 (%.1f%%)", var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", var_expl[2]),
    color = "Group"
  ) +
  theme_bw()




# Heatmap for TOP 50 Species
pc <- 1e-6
X <- mat_f + pc

#Prepare table column
mean_Crohn   <- rowMeans(X[, grp == "Crohn",   drop = FALSE])
mean_Healthy <- rowMeans(X[, grp == "Healthy", drop = FALSE])
log2FC <- log2(mean_Crohn / mean_Healthy)
pval <- apply(X, 1, function(v)
  wilcox.test(v[grp=="Crohn"], v[grp=="Healthy"], exact = FALSE)$p.value
)
padj <- p.adjust(pval, method = "BH")
taxon     <- rownames(mat_f)
direction <- ifelse(log2FC > 0, "Enriched in Crohn", "Depleted in Crohn")
abs_log2FC <- abs(log2FC)

prev_mat <- mat_f > 0
prev_crohn <- rowSums(prev_mat[, grp == "Crohn", drop = FALSE])
prev_ctrl  <- rowSums(prev_mat[, grp == "Healthy", drop = FALSE])

#Build table
res <- tibble(
  taxon = taxon,
  log2FC = log2FC,
  direction = direction,
  pval = pval,
  padj = padj,
  prev_crohn = prev_crohn,
  prev_ctrl = prev_ctrl,
  abs_log2FC = abs_log2FC
)

res <- res[order(res$padj, -res$abs_log2FC), ]

#Select top 50 with presence in 2 samples of at least 1 group
topN <- 50
top_taxa <- res[res$prev_crohn>=2 | res$prev_ctrl>=2,] %>% slice_head(n = topN) %>% pull(taxon)

mat_plot <- log10(mat_f[top_taxa, , drop = FALSE] + pc)

ann_row <- data.frame(Group = grp)
rownames(ann_row) <- colnames(mat_plot)

pheatmap(
  (mat_plot),
  annotation_col = ann_row,
  main = paste0("Top ", topN, " taxa (log10 abundance)\nCrohn vs Healthy; prevalence>=2 in >=1 group"),
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 45,
  border_color = "black"
)


# TOP 25 Enriched/Depleted in Crohn's disease

topN <- 25  
top_crohn<- res[res$prev_crohn>=2,]

enriched_taxa <- top_crohn %>%
  filter(log2FC > 0) %>%                       
  arrange(padj, desc(abs_log2FC), desc(prev_crohn)) %>%
  slice_head(n = topN) %>%
  pull(taxon)

top_healthy<- res[res$prev_ctrl>=2,] 
depleted_taxa <- top_healthy %>%
  filter(log2FC < 0) %>%                       
  arrange(padj, desc(abs_log2FC), desc(prev_ctrl)) %>% 
  slice_head(n = topN) %>%
  pull(taxon)

ann_row <- data.frame(Group = grp)
rownames(ann_row) <- colnames(mat_f)

mat_enriched <- log10(mat_f[enriched_taxa, , drop = FALSE] + pc)
mat_depleted <- log10(mat_f[depleted_taxa, , drop = FALSE] + pc)

pheatmap(
  t(mat_enriched),
  annotation_row = ann_row,
  main = paste0("Top ", topN, " taxa enriched in Crohn's group\n(log10 abundance; prevalence>=2 Crohn)"),
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 315,
  border_color = "black"
)

pheatmap(
  t(mat_depleted),
  annotation_row = ann_row,
  main = paste0("Top ", topN, " taxa depleted in Crohn's group\n(log10 abundance; prevalence>=2 in Healthy)"),
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 315,
  border_color = "black"
)




## Analysis genus Faecalibacterium

mpa_g <- mpa %>%
  filter(str_detect(.data[[tax_col]], "\\|g__[^|]+$")) %>%
  mutate(
    genus = str_extract(.data[[tax_col]], "g__[^|]+$") %>%
      str_remove("^g__") %>%
      str_replace_all("_", " ")
  )

df_fae <- data.frame(
  sample = colnames(mpa_g[,2:6]),
  Faecalibacterium = as.numeric(mpa_g[mpa_g$genus=="Faecalibacterium",2:6]),
  group=grp
)

wt <- wilcox.test(Faecalibacterium ~ group, data = df_fae, exact = FALSE)


p1 <- ggplot(df_fae, aes(x = group, y = Faecalibacterium)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  theme_bw() +
  labs(
    title = "Faecalibacterium abundance (raw)",
    subtitle = paste0("Wilcoxon p = ", signif(wt$p.value, 3))
  )

p2 <- ggplot(df_fae, aes(x = group, y = log10(Faecalibacterium))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, max.overlaps = Inf) +
  theme_bw() +
  labs(
    title = "Faecalibacterium abundance (log10)",
    subtitle = paste0("Wilcoxon p = ", signif(wt$p.value, 3))
  )

print(p1)
print(p2)

