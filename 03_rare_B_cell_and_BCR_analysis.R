# Rare B-cell characterization and BCR analysis: Fig. 2 and Supplementary Fig. 3.
# Includes partial Fig. 2A-D/H code, Fig. 2E-G, and additional DC/repertoire plots.
# Inputs: annotated Seurat objects sobj and b, plus the subtypelist color vector.
# Metadata: Sample, Data, majordisease, Disease, celltype, and subcluster.
# Earlier plotting code is retained; BCR sections include explicit fixes and group filters.
# Create the bcr output directory and load ggpubr for stat_compare_means before running.

library(Seurat)
library(ggplot2)
library(dplyr)
# Fig. 2A: Sample-level Bnaive_IGHV4 proportions.
# The denominator includes all cells in sobj; use a B-cell object for B-cell proportions.
# The current y-axis label names Bnaive_UN and needs correction.
meta <- sobj@meta.data
cell_counts <- sobj@meta.data %>%
  count(Sample,name = 'n')
valid_sample <- cell_counts %>%
  filter(n > 200) %>%
  pull(Sample)
meta <- subset(meta,Sample %in% valid_sample)
meta <- meta[meta$Data != 'GSE242615B',]
qc_counts <- meta %>%
  group_by(sample = meta$Sample, disease = meta$majordisease) %>%
  summarise(n_C1QC = sum(celltype == "Bnaive_IGHV4"),
            total = n(), .groups = "drop") %>%
  mutate(prop_C1QC = n_C1QC / total)
qc_filtered <- subset(qc_counts, disease %in% c('Healthy','KD','MIS-C'))
qc_filtered$disease <- factor(qc_filtered$disease,levels = c('Healthy','KD','MIS-C'))
pdf('bcr/2-2-n.pdf',width = 2.8,height = 5)
ggplot(qc_filtered, aes(x = disease, y = prop_C1QC)) +
  geom_boxplot(color = '#00468B',size = 1,alpha = 0.2,outlier.shape = NA) +
  geom_jitter(color = '#00468B',width = 0.2,size=1.8) +
  theme_bw() +
  ylab("Proportion of Bnaive_UN per sample") +
  xlab("Disease group") +
  theme_classic() +
  stat_compare_means(label = 'p.format',
                     method = 'wilcox.test',
                     ref.group = 'Healthy',
                     bracket.size = 0.2) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size=11),
    axis.text.y = element_text(size=11),
    legend.position = "none",
    axis.ticks.x = element_line(size = 0.9,color = 'black'),
    axis.ticks.y = element_line(size = 0.9,color = 'black')
  ) + coord_cartesian(ylim = c(0, 1.2 * max(qc_filtered$prop_C1QC)))
dev.off()

# Additional analysis: DC proportions (not Fig. 2).
# Plot labels currently refer to Bnaive_4 and do not match the counted DC cells.
meta <- sobj@meta.data
meta <- subset(meta,Sample %in% valid_sample)
qc_counts <- meta %>%
  group_by(sample = meta$Sample, disease = meta$majordisease) %>%
  summarise(n_C1QC = sum(subcluster == "DC"),
            total = n(), .groups = "drop") %>%
  mutate(prop_C1QC = n_C1QC / total)
qc_filtered <- subset(qc_counts, disease %in% c('Healthy','KD','MIS-C'))
qc_filtered$disease <- factor(qc_filtered$disease,levels = c('Healthy','KD','MIS-C'))

ggplot(qc_filtered, aes(x = disease, y = prop_C1QC)) +
  geom_boxplot(color = '#E41A1C',size = 0.6,alpha = 0.2,outlier.shape = NA) +
  geom_jitter(color = '#E41A1C',width = 0.2,size=1) +
  theme_bw() +
  ylab("Proportion of Bnaive_4 per sample") +
  xlab("Disease group") +
  ggtitle("Bnaive_4 proportion per sample in B Cells") +
  theme_classic() +
  stat_compare_means(label = 'p.format',
                     method = 'wilcox.test',
                     ref.group = 'Healthy',
                     bracket.size = 0.25) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size=11),
    axis.text.y = element_text(size=11),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(0, 1.2 * max(qc_filtered$prop_C1QC)))

# Additional analysis: pDC proportions in Healthy, SLE, and DLE groups (not Fig. 2).
meta <- sobj@meta.data
meta <- subset(meta,Sample %in% valid_sample)
qc_counts <- meta %>%
  group_by(sample = meta$Sample, disease = meta$Disease) %>%
  summarise(n_C1QC = sum(celltype == "pDC_LILRA4"),
            total = n(), .groups = "drop") %>%
  mutate(prop_C1QC = n_C1QC / total)
qc_filtered <- subset(qc_counts, disease %in% c('Healthy','SLE','DLE'))
qc_filtered$disease <- factor(qc_filtered$disease,levels = c('Healthy','SLE','DLE'))

ggplot(qc_filtered, aes(x = disease, y = prop_C1QC)) +
  geom_boxplot(color = '#B7EA63',size = 0.6,alpha = 0.2,outlier.shape = NA) +
  geom_jitter(color = '#B7EA63',width = 0.2,size=1) +
  theme_bw() +
  ylab("Proportion of pDC per sample") +
  xlab("Disease group") +
  ggtitle(' ') +
  theme_classic() +
  stat_compare_means(label = 'p.format',
                     method = 'wilcox.test',
                     ref.group = 'Healthy',
                     bracket.size = 0.25) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,size=11),
    axis.text.y = element_text(size=11),
    legend.position = "none"
  )
# Fig. 2B: IGHV expression on the B-cell UMAP.
# Additional plots show TMEM131L expression and B-cell annotations (Supplementary Fig. 3A).
library(Seurat)
library(ggplot2)
library(patchwork)
clean_cols <- c("#e9e5d9", "#6b443f")
minimal_theme <- theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold.italic", color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

p1 <- FeaturePlot(b, features = "IGHV4-31", cols = clean_cols, pt.size = 0.5, order = F) + minimal_theme
p2 <- FeaturePlot(b, features = "IGHV4-28", cols = clean_cols, pt.size = 0.5, order = F) + minimal_theme
p3 <- FeaturePlot(b, features = "IGHV4-61", cols = clean_cols, pt.size = 0.5, order = F) + minimal_theme
p4 <- FeaturePlot(b, features = "IGHV3-20", cols = clean_cols, pt.size = 0.5, order = F) + minimal_theme
p5 <- FeaturePlot(b, features = "IGHV3-43", cols = clean_cols, pt.size = 0.5, order = F) + minimal_theme
p6 <- FeaturePlot(b, features = "TMEM131L", cols = clean_cols, pt.size = 0.5, order = F) + minimal_theme

pdf('bcr/2-f6.pdf',width = 3.5,height = 3)
p6
dev.off()

pdf('bcr/b_dimplot.pdf',width = 6,height = 3.7)
DimPlot(b, group.by = 'celltype', cols = subtypelist)
dev.off()

# Fig. 2C: Marker expression across naive B-cell subsets.
bnaive_markers <- c('TCL1A','IGHD','IL4R','CD38','CD27','XBP1',
                    'IGHV3-20','IGHV3-43','IGHV4-31','IGHV4-28','IGHV4-61',
                    'TMEM131L','PRAG1','NR4A2','FOS','ISG15','IFIT3')
naive_cells <- c('Bnaive_TCL1A','Bnaive_NR4A2','Bnaive_IGHV3','Bnaive_IGHV4','Bnaive_UN','Bnaive_ISG15')
b_naive_subset <- subset(b,celltype %in% naive_cells)
Idents(b_naive_subset) <- "celltype"

Idents(b_naive_subset) <- factor(Idents(b_naive_subset),
                                levels = c('Bnaive_TCL1A','Bnaive_IGHV3',
                                           'Bnaive_IGHV4','Bnaive_UN',
                                           'Bnaive_NR4A2','Bnaive_ISG15'))
pdf('bcr/naive_dot.pdf',width = 10,height = 4.8)
DotPlot(b_naive_subset,
        features = bnaive_markers,
        cols = c("lightgrey", "#b11e24"),
        dot.scale = 7) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled), shape = 21, stroke = 0.5, color = "black") +
  scale_fill_gradient(low = "lightgrey",high = "#b11e24") +
  RotatedAxis() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(title = "Expression Profile of Rare Naive B Subsets",
       x = "Selected Marker Genes",
       y = "Bnaive Subsets")
dev.off()

# Fig. 2D: IGHV4-28 expression across the three rare B-cell subsets.
b_subset <- subset(sobj,celltype %in% c('Bnaive_IGHV3','Bnaive_IGHV4','Bnaive_UN'))
b_subset$celltype <- factor(b_subset$celltype,levels = c('Bnaive_IGHV3','Bnaive_IGHV4','Bnaive_UN'))
plot_data <- FetchData(b_subset, vars = c("IGHV4-28", "celltype"))
pdf('bcr/vl4.pdf',width = 2.8,height = 3.3)
ggplot(plot_data, aes(x = celltype, y = `IGHV4-28`, fill = celltype)) +
  geom_violin(trim = T, color = NA, scale = "width") +
  scale_fill_manual(values = c("Bnaive_IGHV3" = "darkred",
                               "Bnaive_UN" = "#52616B",
                               "Bnaive_IGHV4" = '#00468B'))+
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(size = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45,hjust = 1),
    text = element_text(size = 12)
  ) +
  labs(title = "  ", y = "Expression Level")
dev.off()

# Fig. 2H: Rare B-cell enrichment across MIS-C clinical groups and KD.
# This section reuses meta from earlier code; verify its cell population and disease labels.
library(Seurat)
library(dplyr)
library(tidyverse)
library(reshape2)

# Calculate observed-to-expected ratios from a cell-by-group contingency table.
ROIE <- function(crosstab){
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

# Divide two matrices element by element.
divMatrix <- function(m1, m2){
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}
df <- as.data.frame(table(sobj$Sample))
valid_samples <- as.character(df[df$Freq > 200,]$Var1)
meta <- meta[meta$Sample %in% valid_samples,]
crosstab <- table(meta$celltype,meta$Disease)

roie_mat <- ROIE(crosstab)
roie_long <- reshape2::melt(roie_mat)
colnames(roie_long) <- c("CellType", "Disease", "ROE")
roie_long$ROE_capped <- pmin(roie_long$ROE, 4)
roie_long$ROE_capped <- pmax(roie_long$ROE_capped, 0)

roie_subset <- roie_long %>%
  dplyr::filter(
    CellType %in% c('Bnaive_IGHV3','Bnaive_IGHV4','Bnaive_UN'),
    Disease %in% c("Healthy",'MIS-C Moderate',"MIS-C Severe(SARS-CoV2)","MIS-C Recovered",'KD')
  )
roie_subset$Disease <- factor(roie_subset$Disease,levels = c("Healthy",'MIS-C Moderate',"MIS-C Severe(SARS-CoV2)","MIS-C Recovered",'KD'))
roie_subset$CellType <- factor(roie_subset$CellType,levels = c('Bnaive_UN',"Bnaive_IGHV4",'Bnaive_IGHV3'))
p <- ggplot(roie_subset, aes(x = Disease, y = CellType, fill = ROE_capped)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "lightblue", mid = "white", high = "#bb2b2d",
    midpoint = 1, name = "RO/E"
  ) +
  theme_minimal() +
  labs(title = "", x = "Disease", y = "Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size=10),
    panel.grid = element_blank()
  )
pdf('bcr/roe.pdf',width =5.4 ,height = 4)
p
dev.off()

# Supplementary Fig. 3C: Rare B-cell enrichment across diseases.
# Calculate observed-to-expected ratios from a cell-by-group contingency table.
ROIE <- function(crosstab){
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

# Divide two matrices element by element.
divMatrix <- function(m1, m2){
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

meta <- sobj@meta.data
meta <- meta[meta$majordisease != 'Febrile control+URI' & meta$majordisease != 'idiopathic intracranial hypertension',]
df <- as.data.frame(table(sobj$Sample))
valid_samples <- as.character(df[df$Freq > 200,]$Var1)
meta <- meta[meta$Sample %in% valid_samples,]
crosstab <- table(meta$celltype,meta$Disease)

roie_mat <- ROIE(crosstab)
roie_long <- reshape2::melt(roie_mat)
colnames(roie_long) <- c("CellType", "Disease", "ROE")
roie_long$ROE_capped <- pmin(roie_long$ROE, 3)
roie_long$ROE_capped <- pmax(roie_long$ROE_capped, 0)

roie_long$CellType <- factor(
  roie_long$CellType,
  levels = roie_long %>%
    group_by(CellType) %>%
    summarise(max_roe = max(ROE, na.rm = TRUE)) %>%
    arrange(desc(max_roe)) %>%
    pull(CellType)
)

roie_long$Disease <- factor(
  roie_long$Disease,
  levels = roie_long %>%
    group_by(Disease) %>%
    summarise(mean_roe = mean(ROE, na.rm = TRUE)) %>%
    arrange(desc(mean_roe)) %>%
    pull(Disease)
)
roie_subset <- roie_long %>%
  dplyr::filter(
    CellType %in% c("Bnaive_UN", "Bnaive_IGHV3",'Bnaive_IGHV4'),
  )
roie_subset$Disease <- factor(roie_subset$Disease,levels = c('MIS-C',"KD",'Healthy',"Systemic lupus erythematosus",
                                                             "Behçet’s disease",'CD','primary Sjögren’s syndrome','Co-morbid CD and Spondyloarthritis',
                                                             'Giant cell arteritis ','Mutiple sclerosis','RA','SLE','Spondyloarthritis'))
roie_subset$CellType <- factor(roie_subset$CellType,levels = c("Bnaive_UN", "Bnaive_IGHV3",'Bnaive_IGHV4'))
pdf('bcr/3rare_roe_add.pdf',width = 6,height = 3.5)
ggplot(roie_subset, aes(x = Disease, y = CellType, fill = ROE_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", ROE_capped)), size = 0) +
  scale_fill_gradient2(
    low = "lightblue", mid = "white", high = "#bb2b2d",
    midpoint = 1., name = "RO/E"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal() +
  labs(x = "Disease", y = "Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(0, 0, 0, 0)
  )
dev.off()

# Re-export the same cross-disease heatmap; this block does not recompute clinical-group enrichment.
pdf('bcr/3rare_roe_progression.pdf',width = 5,height = 3.5)
ggplot(roie_subset, aes(x = Disease, y = CellType, fill = ROE_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", ROE_capped)), size = 0) +
  scale_fill_gradient2(
    low = "lightblue", mid = "white", high = "#bb2b2d",
    midpoint = 1., name = "RO/E"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal() +
  labs(x = "Disease", y = "Cell Type") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    plot.margin = margin(0, 0, 0, 0)
  )
dev.off()

# BCR analysis: Fig. 2E-G and additional repertoire plots.
# Input: sobj with CTgene/CTnt annotations and a precomputed bcr_downsampled table.
# V-J pairings are gene-usage categories, not individual clonotypes.
library(scRepertoire)
library(tidyr)
library(stringr)

bcr_data <- sobj@meta.data %>%
  filter(celltype %in% c("Bnaive_IGHV3", "Bnaive_IGHV4", "Bnaive_UN"))

# Parse V and J genes using the original CTgene format.
bcr_processed <- bcr_data %>%
  filter(!is.na(CTgene) & CTgene != "None") %>%
  separate(CTgene, c("v_gene", "temp_rest"), sep = "\\.",
           extra = "merge", remove = FALSE) %>%
  mutate(j_gene = str_extract(CTgene, "IGHJ[0-9]+")) %>%
  select(-temp_rest)
head(bcr_processed[, c("CTgene", "v_gene", "j_gene")])

# The original downsampling step was not supplied; do not substitute pooled data.
if (!exists("bcr_downsampled")) {
  stop("Provide bcr_downsampled from the original downsampling analysis before running Fig. 2E-G.")
}
stopifnot(all(c("majordisease", "v_gene", "j_gene") %in% names(bcr_downsampled)))

# Additional plot: V-J pairing frequencies by disease.
vj_summary <- bcr_downsampled %>%
  group_by(majordisease, v_gene, j_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(majordisease) %>%
  mutate(proportion = count / sum(count))
pdf("BCR_VJ_frequency_heatmap.pdf", width = 8, height = 4)
ggplot(vj_summary, aes(x = v_gene, y = j_gene, fill = proportion)) +
  geom_tile() + facet_wrap(~majordisease) +
  scale_fill_viridis_c(option = "magma") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()

# Rank V-J pairings by frequency difference; no significance test is performed.
comparison <- bcr_downsampled %>%
  filter(!is.na(v_gene) & v_gene != "None", !is.na(j_gene) & j_gene != "None") %>%
  group_by(majordisease, v_gene, j_gene) %>%
  tally() %>%
  group_by(majordisease) %>%
  mutate(freq = n / sum(n))
diff_df <- comparison %>%
  select(-n) %>%
  pivot_wider(names_from = majordisease, values_from = freq, values_fill = 0) %>%
  mutate(diff = `MIS-C` - Healthy) %>%
  arrange(desc(diff))
head(diff_df, 20)

# Additional plot: the 15 largest frequency differences.
top_diff <- head(diff_df, 15) %>%
  unite("vj", v_gene, j_gene, remove = FALSE)
pdf("BCR_VJ_frequency_difference_barplot.pdf")
ggplot(top_diff, aes(x = reorder(vj, diff), y = diff, fill = diff > 0)) +
  geom_bar(stat = "identity") + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  theme_classic() +
  labs(title = "V-J Pairings Enriched in MIS-C", x = "V-J Combination",
       y = "Difference in Proportion (MIS-C - Healthy)") +
  theme(legend.position = "none")
dev.off()

# Fig. 2E: frequency differences for the top 20 V-J pairings.
# Expansion Index = MIS-C frequency minus Healthy frequency.
# Outlines indicate absence in the analyzed Healthy data, not biological exclusivity.
plot_df <- diff_df %>%
  head(20) %>%
  unite("vj", v_gene, j_gene, remove = FALSE) %>%
  mutate(is_new = ifelse(Healthy == 0, "New in Disease", "Enriched"))
pdf("hihiBCR_VJ_frequency_difference_barplot.pdf")
ggplot(plot_df, aes(x = reorder(vj, diff), y = diff)) +
  geom_segment(aes(xend = reorder(vj, diff), y = 0, yend = diff),
               color = "#E0E6ED", size = 0.8) +
  geom_point(aes(color = diff, size = `MIS-C`), alpha = 0.9) +
  geom_point(data = filter(plot_df, is_new == "New in Disease"),
             aes(x = vj, y = diff), color = "#A020F0", size = 5, shape = 1, stroke = 1.2) +
  scale_color_gradientn(colors = c("#2B2D42", "#3A86FF", "#8338EC", "#00F5FF")) +
  scale_size_continuous(range = c(3, 8), breaks = c(0.05, 0.1, 0.2)) +
  coord_flip() + theme_minimal() +
  labs(title = "V-J Pairing Frequencies in MIS-C", x = NULL,
       y = "Frequency Delta (MIS-C vs Healthy)", size = "Abundance in MIS-C",
       color = "Expansion Index") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "#F0F0F0", linetype = "dashed"),
        plot.title = element_text(face = "bold", size = 15, color = "#1A1A1A"),
        plot.subtitle = element_text(size = 10, color = "#666666", face = "italic"),
        axis.text = element_text(color = "#333333", size = 10, family = "sans"),
        axis.title.x = element_text(size = 11, margin = margin(t = 10)),
        legend.position = "right", legend.background = element_rect(fill = "white", color = NA))
dev.off()

# Fig. 2F: subset distribution of the top 15 pairings in the original BCR data.
diff_df <- diff_df %>% mutate(vj = paste0(v_gene, "_", j_gene))
top_vj_list <- head(diff_df$vj, 15)
subgroup_analysis <- bcr_processed %>%
  unite("vj", v_gene, j_gene, remove = FALSE) %>%
  filter(vj %in% top_vj_list, majordisease %in% c("Healthy", "MIS-C")) %>%
  group_by(celltype, majordisease, vj) %>%
  summarise(cell_count = n(), .groups = "drop")
pdf("Fig2F_VJ_subset_distribution.pdf", width = 4.5, height = 5.5)
ggplot(subgroup_analysis, aes(x = celltype, y = reorder(vj, cell_count), size = cell_count)) +
  geom_point(aes(color = majordisease), alpha = 0.8) +
  scale_color_manual(values = c("Healthy" = "#95A5A6", "MIS-C" = "#00D2FF")) +
  scale_size_continuous(range = c(2, 10)) + facet_wrap(~majordisease) + theme_minimal() +
  labs(x = "Bnaive Subsets", y = "Top Enriching V-J Pairs", size = "Cell Count", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "#FAFAFA", color = NA),
        strip.text = element_text(face = "bold", size = 12))
dev.off()

# Fig. 2G: within-subset V-J frequencies in MIS-C only.
# The supplied code pooled diseases; this filter matches the figure legend.
dist_data <- bcr_processed %>%
  filter(majordisease == "MIS-C") %>%
  unite("vj", v_gene, j_gene, remove = FALSE) %>%
  group_by(celltype, vj) %>% tally() %>%
  mutate(freq = n / sum(n)) %>% ungroup()
top_10_vj <- dist_data %>%
  group_by(vj) %>% summarise(total_n = sum(n), .groups = "drop") %>%
  top_n(10, total_n) %>% pull(vj)
dist_data <- dist_data %>% mutate(vj_group = ifelse(vj %in% top_10_vj, vj, "Others"))
pdf("Fig2G_MISC_VJ_composition.pdf", width = 4, height = 6)
ggplot(dist_data, aes(x = celltype, y = freq, fill = vj_group)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", size = 0.1) +
  scale_fill_manual(values = c("Others" = "#ECF0F1",
    setNames(scales::hue_pal(l = 60)(length(top_10_vj)), top_10_vj))) +
  theme_minimal() +
  labs(title = "V-J Pairing Composition in MIS-C", x = NULL,
       y = "Proportion of Repertoire", fill = "V-J Combination") +
  theme(panel.grid = element_blank(), axis.line.y = element_line(color = "black"),
        axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1),
        legend.position = "right", panel.border = element_rect(color = "black", fill = NA, size = 1))
dev.off()

# Additional analysis: CTnt frequency classes pooled within each subset and disease.
# Counts are not donor-specific and should not be interpreted as within-donor expansion.
library(patchwork)
bcr_summary <- bcr_data %>%
  filter(!is.na(CTnt) & CTnt != "") %>%
  group_by(celltype, majordisease, CTnt) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(cloneSize = case_when(n == 1 ~ "Singleton (1)", n == 2 ~ "Small (2)", n > 2 ~ "Expanded (>2)"))
plot_data <- bcr_summary %>%
  group_by(celltype, majordisease, cloneSize) %>%
  summarise(count = sum(n), .groups = "drop") %>%
  group_by(celltype, majordisease) %>% mutate(proportion = count / sum(count))
pdf("BCR_pooled_CTnt_frequency_classes.pdf", width = 3, height = 4)
ggplot(plot_data, aes(x = majordisease, y = proportion, fill = cloneSize)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") + facet_wrap(~celltype) +
  scale_fill_manual(values = c("Singleton (1)" = "#00D2FF", "Small (2)" = "#0099CC", "Expanded (>2)" = "#34495E")) +
  theme_bw() + labs(title = "Pooled CTnt Frequency Classes", x = "Disease Group",
                   y = "Proportion of B cells", fill = "Clonal Size") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
