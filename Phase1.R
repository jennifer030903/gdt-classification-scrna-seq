library(Seurat)
library(dplyr)
library(tidyr)
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(pROC)
library(caret)
library(rsample)
library(tibble)
library(PRROC)
library(Metrics)


# ----------------------------------------------------------
# Call data
# ----------------------------------------------------------

Sys.setenv(LANG = "en")
seurat_obj <- readRDS("C:/Users/jysje/OneDrive/Desktop/Uni/4-1/BBMS4001/gdtcells.rds")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)


# make a vector for the six genes 
six_genes <- c("TRDC", "TRGC1", "TRGC2", "TRAC", "TRBC1", "TRBC2")

# ----------------------------------------------------------
# Double checking that the scale data exists 
# ----------------------------------------------------------

# Checking a few values 
seurat_obj@assays$RNA@scale.data[1:10, 1:10]


# check if all six genes are in the scale data 
six_genes %in% rownames(seurat_obj@assays$RNA@scale.data)

# the mean and sd aren't scaled well 
for (gene in six_genes) {
  gene_mean <- mean(seurat_obj@assays$RNA@scale.data[gene, ])
  gene_sd   <- sd(seurat_obj@assays$RNA@scale.data[gene, ])
  
  if (abs(gene_mean - 0) < 1e-6) {
    print(paste0("mean for ", gene, " ≈ 0"))
  } else {
    print(paste0("mean for ", gene, " = ", round(gene_mean, 3)))
  }
  
  if (abs(gene_sd - 1) < 1e-6) {
    print(paste0("sd for ", gene, " ≈ 1"))
  } else {
    print(paste0("sd for ", gene, " = ", round(gene_sd, 3)))
  }
}

# ----------------------------------------------------------
# Scaling the data with all genes 
# ----------------------------------------------------------

all.genes <- rownames(seurat_obj)

# Additional line with assay = "RNA" to specify where the data should be stored
seurat_obj <- ScaleData(seurat_obj, features = all.genes, assay = "RNA")

# ----------------------------------------------------------
# Perform Linear Dimensional Reduction
# ----------------------------------------------------------

# Each cell gets a score for each PC, stored in: seurat_obj@reductions$pca@cell.embeddings
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# To check how many PCs there are 
dim(seurat_obj@reductions$pca@cell.embeddings)

# Get the PCA embeddings (cell scores) top 20 
head(seurat_obj@reductions$pca@cell.embeddings[, 1:20])

# Get the total number of PCs 
ncol(seurat_obj@reductions$pca@cell.embeddings)

# Get the total number of cells 
nrow(seurat_obj@reductions$pca@cell.embeddings)

print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
DimHeatmap(seurat_obj, dims = 1: 18, cells = 500, balanced = TRUE)
ElbowPlot(seurat_obj)

# ----------------------------------------------------------
# Cluster the cells & UMAP
# ----------------------------------------------------------

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:18, assay = "RNA")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

#view the names and which ones are stored
names(seurat_obj@graphs)

# See the actual graph 
seurat_obj@graphs$RNA_nn[1:5, 1:5]
seurat_obj@graphs$RNA_snn[1:5, 1:5]


# Making the UMAP 
seurat_obj <- RunUMAP(seurat_obj, dims = 1:18)
DimPlot(seurat_obj, reduction = "umap")


# ----------------------------------------------------------
# Finding Differentially expressed features (cluster biomarkers)
# ----------------------------------------------------------

seurat.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
seurat.markers |> 
  group_by(cluster) |> 
  filter(avg_log2FC > 1)

seurat_obj@misc$markers <- seurat.markers
saveRDS(suerat_obj, file = "gtcells_newMarkerData.rds")


# Plot violin plot 
VlnPlot(seurat_obj, features = six_genes)

# plot the raw data 
VlnPlot(seurat_obj, features = six_genes, slot = "counts", log = TRUE)

# Feature plot showing UMAP of the 6 genes 
FeaturePlot(seurat_obj, features = six_genes)

# DoHeatmap() 
# This data is stored in pbmc@misc$top10_markers 
seurat_obj.markers |> 
  group_by(cluster) |> 
  filter(avg_log2FC > 1) |> 
  slice_head(n = 10) |> 
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()

# DimPlot() 
# Extract the predicted labels from metadata
new.cluster.ids <- seurat_obj@meta.data$predicted.celltype.l2

# Assign names to match cluster levels 
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

seurat_obj@meta.data$seurat_clusters