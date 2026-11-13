# Fig. 3A-F and Supplementary Fig. 4A-F.
# Input: sobj containing all annotated CD4 T cells, normalized RNA data, and UMAP.
# Required metadata: Sample, majordisease, celltype.
# Abundance analyses retain samples with >200 CD4 T cells.
# Expression and UCell analyses use cells directly, as in the source scripts.

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(msigdbr)
library(UCell)
library(ComplexHeatmap)

out_dir <- "fig3_treg"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)
cd4 <- sobj
DefaultAssay(cd4) <- "RNA"
stopifnot(all(c("Sample", "majordisease", "celltype") %in% colnames(cd4[[]])))

# Standardize legacy subset and disease labels.
cell_map <- c(CD4Treg_Naive = "CD4rTreg", rTreg = "CD4rTreg",
              CD4Treg_Act = "CD4rTreg-Primed", rTreg_Primed = "CD4rTreg-Primed",
              CD4Treg_Eff = "CD4aTreg", aTreg = "CD4aTreg")
disease_map <- c("Behçet’s disease" = "BD", "Behçet's disease" = "BD",
                 "Systemic lupus erythematosus" = "SLE",
                 "idiopathic intracranial hypertension" = "IIH",
                 "Giant cell arteritis" = "GCA", "Mutiple sclerosis" = "MS",
                 "Multiple sclerosis" = "MS", "Febrile control+URI" = "FEB",
                 "Co-morbid CD and Spondyloarthritis" = "CD-axSpA",
                 CDAS = "CD-axSpA", AS = "axSpA")
rename_labels <- function(x, map) {
  x <- trimws(as.character(x))
  hit <- x %in% names(map)
  x[hit] <- unname(map[x[hit]])
  x
}
cd4$celltype <- rename_labels(cd4$celltype, cell_map)
cd4$majordisease <- rename_labels(cd4$majordisease, disease_map)
treg_types <- c("CD4rTreg", "CD4rTreg-Primed", "CD4aTreg")
cd4$subcluster <- ifelse(cd4$celltype %in% treg_types, "CD4Treg", "Non-Treg")
meta_all <- cd4[[]]
stopifnot(!anyNA(meta_all[, c("Sample", "majordisease", "celltype")]))
valid_samples <- meta_all %>% count(Sample) %>% filter(n > 200) %>% pull(Sample)
meta <- meta_all %>% filter(Sample %in% valid_samples)
stopifnot(nrow(meta) > 0)

save_plot <- function(p, name, width = 5, height = 4) {
  ggsave(file.path(out_dir, name), plot = p, width = width, height = height)
}
box_plot <- function(data, value, label, color = "#8aace4") {
  ggplot(data, aes(x = majordisease, y = .data[[value]])) +
    geom_boxplot(color = color, outlier.shape = NA, width = 0.6) +
    geom_jitter(color = color, width = 0.2, size = 1) +
    theme_classic() + labs(x = NULL, y = label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Fig. 3A: all three Treg subsets divided by total CD4 T cells per sample.
abundance <- meta %>% group_by(Sample, majordisease) %>%
  summarise(total_CD4 = n(), n_Treg = sum(subcluster == "CD4Treg"),
            n_aTreg = sum(celltype == "CD4aTreg"),
            n_Th17 = sum(celltype == "CD4Tem_Th17"), .groups = "drop") %>%
  mutate(Treg_prop = n_Treg / total_CD4, aTreg_prop = n_aTreg / total_CD4,
         Th17_prop = n_Th17 / total_CD4)
disease_order <- abundance %>% group_by(majordisease) %>%
  summarise(median_prop = median(Treg_prop), .groups = "drop") %>%
  arrange(median_prop) %>% pull(majordisease)
abundance$majordisease <- factor(abundance$majordisease, levels = disease_order)
p <- box_plot(abundance, "Treg_prop", "Proportion of Tregs among CD4 T cells") +
  stat_compare_means(method = "kruskal.test", label = "p.format")
save_plot(p, "Fig3A_Treg_proportions.pdf", 6, 4)
write.csv(abundance, file.path(out_dir, "sample_cell_proportions.csv"), row.names = FALSE)

# Fig. 3B: marker expression in the three Treg subsets.
tregs <- subset(cd4, cells = rownames(meta_all)[meta_all$subcluster == "CD4Treg"])
tregs$celltype <- factor(tregs$celltype, levels = treg_types)
Idents(tregs) <- "celltype"
markers <- c("SELL", "CCR7", "TCF7", "LEF1", "CXCR4", "CXCR3", "ITGA4",
             "S1PR1", "FOXP3", "IL2RA", "CTLA4", "TIGIT")
p <- DotPlot(tregs, features = markers, group.by = "celltype", cols = "RdYlBu") +
  coord_flip() + theme_classic() + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "Fig3B_Treg_marker_expression.pdf", 3.5, 4)

# Fig. 3C: pooled cell proportions within the Treg compartment by disease.
selected_diseases <- c("Healthy", "CD", "axSpA", "KD", "RA", "SLE", "BD", "pSS")
composition <- meta %>%
  filter(subcluster == "CD4Treg", majordisease %in% selected_diseases) %>%
  count(majordisease, celltype) %>%
  complete(majordisease, celltype = treg_types, fill = list(n = 0)) %>%
  group_by(majordisease) %>% mutate(percent = n / sum(n) * 100) %>% ungroup()
composition_order <- composition %>% filter(celltype == "CD4aTreg") %>%
  arrange(percent) %>% pull(majordisease)
composition$majordisease <- factor(composition$majordisease, levels = composition_order)
composition$celltype <- factor(composition$celltype, levels = treg_types)
treg_colors <- c(CD4rTreg = "#5e6683", "CD4rTreg-Primed" = "#f4dfd8", CD4aTreg = "#97c9ab")
p <- ggplot(composition, aes(majordisease, percent, fill = celltype)) +
  geom_col(width = 0.7) + scale_fill_manual(values = treg_colors) +
  scale_y_continuous(expand = c(0, 0)) + theme_classic() +
  labs(x = NULL, y = "Percentage (%)", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "Fig3C_Treg_composition.pdf", 5, 4)

# Fig. 3D: disease-level mean UCell scores, standardized across diseases.
# Retain the original gene sets; verify the resting-versus-activated signature direction.
get_genes <- function(category, name, subcategory = NULL) {
  tab <- if (is.null(subcategory)) msigdbr(species = "Homo sapiens", category = category) else
    msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory)
  genes <- unique(tab$gene_symbol[tab$gs_name == name])
  if (!length(genes)) stop(paste("Gene set not found:", name))
  genes
}
signatures <- list(
  Resting_vs_Activated_UP = get_genes("C7", "GSE15659_RESTING_VS_ACTIVATED_TREG_UP"),
  IL2_STAT5 = get_genes("H", "HALLMARK_IL2_STAT5_SIGNALING"),
  TGF_beta = get_genes("H", "HALLMARK_TGF_BETA_SIGNALING"),
  OXPHOS = get_genes("H", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  Migration = get_genes("C5", "GOBP_LEUKOCYTE_MIGRATION", "BP"))
atregs <- subset(tregs, celltype == "CD4aTreg")
atregs <- AddModuleScore_UCell(atregs, features = signatures)
score_cols <- paste0(names(signatures), "_UCell")
score_means <- atregs[[]] %>% group_by(majordisease) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop") %>%
  tibble::column_to_rownames("majordisease")
colnames(score_means) <- c("Resting vs activated Treg UP", "IL-2/STAT5 signaling",
                          "TGF-beta signaling", "Oxidative phosphorylation", "Cell migration")
plot_matrix <- t(scale(score_means))
heatmap_order <- c(intersect("Healthy", colnames(plot_matrix)),
                   sort(setdiff(colnames(plot_matrix), "Healthy")))
plot_matrix <- plot_matrix[, heatmap_order, drop = FALSE]
pdf(file.path(out_dir, "Fig3D_aTreg_signature_scores.pdf"), width = 7, height = 3.7)
draw(Heatmap(plot_matrix, name = "Z-score",
             col = colorRampPalette(c("#8fb4be", "white", "#da9599"))(100),
             cluster_columns = FALSE, cluster_rows = TRUE,
             rect_gp = grid::gpar(col = "white", lwd = 2), column_names_rot = 45))
dev.off()

# Fig. 3E: cell-level Th17 signature scores in aTregs, Healthy versus pSS.
th17_genes <- get_genes("C7", "GSE14308_TH17_VS_NATURAL_TREG_UP")
atregs <- AddModuleScore_UCell(atregs, features = list(Th17 = th17_genes))
score_data <- atregs[[]] %>% filter(majordisease %in% c("Healthy", "pSS")) %>%
  mutate(majordisease = factor(majordisease, levels = c("Healthy", "pSS")))
group_colors <- c(Healthy = "#999999", pSS = "#ea9e58")
p <- ggplot(score_data, aes(majordisease, Th17_UCell, fill = majordisease, color = majordisease)) +
  geom_violin(alpha = 0.2, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, color = "black", fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = group_colors) + scale_color_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() + labs(x = NULL, y = "Th17 UCell score in aTregs") +
  theme(legend.position = "none")
save_plot(p, "Fig3E_Th17_signature_in_aTregs.pdf", 3, 3.5)

# Fig. 3F: sample-level Pearson correlations with linear regression fits.
cor_data <- abundance %>% mutate(Th17_percent = Th17_prop * 100, aTreg_percent = aTreg_prop * 100)
for (group in c("Healthy", "pSS")) {
  dat <- cor_data %>% filter(majordisease == group)
  p <- ggscatter(dat, x = "Th17_percent", y = "aTreg_percent",
                 color = unname(group_colors[group]), alpha = 0.6, size = 4,
                 add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson") +
    theme_classic() + labs(title = group, x = "Th17 abundance (%)", y = "aTreg abundance (%)")
  save_plot(p, paste0("Fig3F_Th17_aTreg_", group, ".pdf"), 3.5, 3.5)
}

# Supplementary Fig. 4A: total Treg proportions in selected diseases.
dat <- abundance %>% filter(majordisease %in% c("Healthy", "KD", "pSS", "SLE")) %>%
  mutate(majordisease = factor(majordisease, levels = c("Healthy", "KD", "pSS", "SLE")))
p <- box_plot(dat, "Treg_prop", "Proportion of Tregs among CD4 T cells") +
  stat_compare_means(method = "wilcox.test", ref.group = "Healthy", label = "p.format")
save_plot(p, "SFig4A_Treg_proportions.pdf", 3.4, 4)

# Supplementary Fig. 4B-C: FOXP3 expression and CD4 T-cell annotations.
save_plot(FeaturePlot(cd4, features = "FOXP3"), "SFig4B_FOXP3_UMAP.pdf", 4, 4)
save_plot(DimPlot(cd4, group.by = "celltype"), "SFig4C_CD4_subsets.pdf", 7, 4)

# Supplementary Fig. 4D: six paired values supplied in the source script.
# Patient labels are placeholders; the original per-patient extraction was not supplied.
before <- c(0.005617978, 0.014015417, 0.005213270, 0.024291498, 0.024890190, 0.015968064)
after <- c(0.018714402, 0.016724519, 0.034988713, 0.018348624, 0.047752809, 0.018832392)
paired_test <- wilcox.test(before, after, paired = TRUE, alternative = "less")
paired_data <- data.frame(Patient = rep(paste0("P", 1:6), 2),
                          Time = factor(rep(c("Before", "After"), each = 6), levels = c("Before", "After")),
                          Value = c(before, after))
p <- ggplot(paired_data, aes(Time, Value, group = Patient)) +
  geom_line(color = "grey80", linetype = "dashed") +
  geom_point(aes(color = Time), size = 4) +
  scale_color_manual(values = c(Before = "#51999f", After = "#97c9ab")) +
  theme_classic() + labs(x = NULL, y = "Proportion of aTregs among CD4 T cells",
                         subtitle = paste0("One-sided paired Wilcoxon, p = ", signif(paired_test$p.value, 3))) +
  theme(legend.position = "none")
save_plot(p, "SFig4D_paired_IVIG_aTreg.pdf", 4, 4)

# Supplementary Fig. 4E: retain the source script's aTreg-only expression analysis.
# This panel does not show all three Treg subsets pooled together.
expr_cells <- subset(atregs, majordisease %in% selected_diseases)
expr_cells$majordisease <- factor(expr_cells$majordisease,
  levels = c("Healthy", "CD", "axSpA", "BD", "RA", "SLE", "KD", "pSS"))
p <- VlnPlot(expr_cells, features = c("FOXP3", "CTLA4", "CCR4"),
             group.by = "majordisease", pt.size = 0, ncol = 3) & NoLegend()
save_plot(p, "SFig4E_aTreg_marker_expression.pdf", 12, 3.2)

# Supplementary Fig. 4F: Th17 proportions among CD4 T cells.
dat <- abundance %>% filter(majordisease %in% c("Healthy", "pSS")) %>%
  mutate(majordisease = factor(majordisease, levels = c("Healthy", "pSS")))
p <- box_plot(dat, "Th17_prop", "Proportion of Th17 cells among CD4 T cells", "#41b5c4") +
  stat_compare_means(method = "wilcox.test", label = "p.format")
save_plot(p, "SFig4F_Th17_proportions.pdf", 3, 4)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
