library(Seurat)
library(ggplot2)
library(ggrepel)
library(harmony)

sobj[["percent.mt"]] <- PercentageFeatureSet(
  sobj,
  assay = "RNA",
  pattern = "^MT-"
)

# Filter low-quality cells
sobj <- subset(
  sobj,
  subset = nFeature_RNA > 500 & percent.mt < 20
)

# Split RNA layers by patient
sobj[["RNA"]] <- split(sobj[["RNA"]], f = sobj$Patient)

sobj <- NormalizeData(sobj)
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj)
sobj <- RunPCA(sobj)
sobj <- JoinLayers(sobj)

sobj <- FindNeighbors(sobj, dims = 1:20)
sobj <- FindClusters(sobj, resolution = 1.2)
sobj <- RunUMAP(sobj, dims = 1:20)
DimPlot(sobj,label = T)


sobj <- RunHarmony(
  sobj,
  group.by.vars='Patient',
  reduction.use = "pca",
  dims.use = NULL,
  reduction.save = "harmony",
  project.dim=F,
  lambda = 0.8
)
sobj <- FindNeighbors(sobj, reduction = "harmony", dims = 1:20)
sobj <- FindClusters(sobj, resolution = 1.2)
sobj <- RunUMAP(sobj, reduction = "harmony", dims = 1:20)

p <- DimPlot(sobj) + NoLegend()
coord <- Embeddings(sobj) %>% as.data.frame()
coord$celltype <- sobj$celltype
label_coord <- aggregate(coord[,1:2],by= list(coord$celltype),FUN = mean)
colnames(label_coord) <- c('celltype','umap_1','umap_2')


DimPlot(sobj,group.by = 'Data')
DimPlot(sobj,group.by = 'Disease')
DimPlot(sobj,group.by = 'celltype',label = T)

DimPlot(sobj,label = T)
DimPlot(sobj,group.by = 'majorcluster',label = T,cols = cluster_colors) +NoLegend()

library(dplyr)
cluster_colors <- c('B'='#d6bbc4','Plasma'='#ffe4b9','Mono/Macro/Neutro' = '#B2D8CE','DC' = '#d0d8eb','NK'='#e6d6e9','T'='#939eb2')

df <- sobj@meta.data
df <- df[df$Data != "GSE242615B", ]
df_summary <- df %>%
  group_by(Disease, majorcluster) %>%
  tally() %>%
  group_by(Disease) %>%
  mutate(percentage = n / sum(n) * 100)

ggplot(df_summary, aes(x = Disease, y = percentage, fill = majorcluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage (%)", x = "Disease", title = "Composition of Major Clusters by Disease") +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),
        axis.text.y = element_text(size=10))
colors_list <- c('#FDF7F4','#F2ED9A','#B7EA63','#46855E','#2B3455','#8EB486','#997C70','#685752',
                 '#DDE7F2','#DDD4FF','#9488D3','#DFF4F3','#B7E0FF','#FFF5CD','#FFCFB3','#E78F81',
                 '#EBE8DB','#F0F5F9','#C9D6DF','#52616B','#ECF9FF','#F8CBA6','#89BBDF','#878ECD',
                 '#EEEEEE','#D9F9F4','#9CD9DE','#CADAD2','#16A23C','#DCD7C9','#E5E1DA','#C890A7',
                 '#FFC9AD','#F99B36','#86FFD6','#0D747A','#FBF5E5','#9AACCB','#F7EED3','#C9FFF9',
                 '#ABE4EF','#95D4FF','#C3E1DC','#8ABBC1','#FADFA1','#FFE2E2','#FFB9B9','#B8E0EF',
                 '#58B8D1','#EBFFD8','#E3D095','#D5E5D5','#ADB2D4','#EEF1DA','#42A197')

