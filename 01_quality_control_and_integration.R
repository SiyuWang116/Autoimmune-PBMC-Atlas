
#0. Integration
library(Seurat)

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
sobj <- subset(sobj, subset = nFeature_RNA > 500 & percent.mt < 20)

sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$Patient)

sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)

sobj <- IntegrateLayers(
  object = sobj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:50)
sobj <- FindClusters(sobj, resolution = 1.0)
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:50)


#1. UMAP & Feature plot
library(ggplot2)
library(RColorBrewer)

cluster_colors <- c(
  "B" = "#d6bbc4",
  "Plasma" = "#ffe4b9",
  "Mono/Macro/Neutro" = "#B2D8CE",
  "DC" = "#d0d8eb",
  "NK" = "#e6d6e9",
  "T" = "#939eb2"
)

colors_list <- c(
  '#FDF7F4','#F2ED9A','#B7EA63','#46855E','#2B3455',
  '#8EB486','#997C70','#685752','#DDE7F2','#DDD4FF',
  '#9488D3','#DFF4F3','#B7E0FF','#FFF5CD','#FFCFB3',
  '#E78F81','#EBE8DB','#F0F5F9','#C9D6DF','#52616B',
  '#ECF9FF','#F8CBA6','#89BBDF','#878ECD','#EEEEEE',
  '#D9F9F4','#9CD9DE','#CADAD2','#16A23C','#DCD7C9',
  '#E5E1DA','#C890A7','#FFC9AD','#F99B36','#86FFD6',
  '#0D747A','#FBF5E5','#9AACCB','#F7EED3','#C9FFF9',
  '#ABE4EF','#95D4FF','#C3E1DC','#8ABBC1','#FADFA1',
  '#FFE2E2','#FFB9B9','#B8E0EF','#58B8D1','#EBFFD8',
  '#E3D095','#D5E5D5','#ADB2D4','#EEF1DA','#42A197'
)

colors_list <- grDevices::colorRampPalette(colors_list)(
  length(unique(sobj$celltype))
)

DimPlot(sobj, group.by = "Data")
DimPlot(sobj, group.by = "Disease")
DimPlot(sobj, group.by = "celltype", label = TRUE, cols = colors_list)
DimPlot(sobj, group.by = "majorcluster", label = TRUE,
        cols = cluster_colors) + NoLegend()

FeaturePlot(sobj, features = c(
  "CD3D", "CD4", "CD8A", "MS4A1", "CD79A",
  "NKG7", "KLRF1", "LYZ", "FCGR3A", "ITGAX"
))


#2. Cell composition
library(dplyr)

df <- sobj@meta.data %>%
  filter(Data != "GSE242615B") %>%
  count(Disease, majorcluster) %>%
  group_by(Disease) %>%
  mutate(percentage = n / sum(n) * 100)

ggplot(df, aes(x = Disease, y = percentage, fill = majorcluster)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cluster_colors) +
  labs(x = NULL, y = "Percentage (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#3. Cell composition by sample
df <- sobj@meta.data %>%
  filter(Data != "GSE242615B") %>%
  count(Sample, Disease, majorcluster) %>%
  group_by(Sample) %>%
  mutate(percentage = n / sum(n) * 100)

ggplot(df, aes(x = Disease, y = percentage, color = majorcluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.3) +
  scale_color_manual(values = cluster_colors) +
  theme_classic() +
  facet_grid(. ~ majorcluster, scales = "free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())
