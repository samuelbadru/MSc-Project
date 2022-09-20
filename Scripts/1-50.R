##Setup the Seurat Object

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
onefifty.data <- Read10X(data.dir = "C:/Users/User/Desktop/Cellranger Outputs/count/1_50GEXFB/filtered_feature_bc_matrix1_50/Zipped/")

# Initialize the Seurat object with the raw (non-normalized data).
onefiftyGEX <- CreateSeuratObject(counts = onefifty.data[["Gene Expression"]], project = "1-50")
onefiftyGEX

##Standard pre-processing workflow

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
onefiftyGEX[["percent.mt"]] <- PercentageFeatureSet(onefiftyGEX, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(onefiftyGEX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(onefiftyGEX, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(onefiftyGEX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

onefiftyGEX <- subset(onefiftyGEX, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##Normalising the data

onefiftyGEX <- NormalizeData(onefiftyGEX, normalization.method = "LogNormalize", scale.factor = 10000)

##Identification of highly variable features (feature selection)

onefiftyGEX <- FindVariableFeatures(onefiftyGEX, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(onefiftyGEX), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(onefiftyGEX)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

##Scaling the data

all.genes <- rownames(onefiftyGEX)
onefiftyGEX <- ScaleData(onefiftyGEX, features = all.genes)

##Perform linear dimensional reduction

onefiftyGEX <- RunPCA(onefiftyGEX, features = VariableFeatures(object = onefiftyGEX))

# Examine and visualize PCA results a few different ways
print(onefiftyGEX[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(onefiftyGEX, dims = 1:2, reduction = "pca")

DimPlot(onefiftyGEX, reduction = "pca")

DimHeatmap(onefiftyGEX, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(onefiftyGEX, dims = 1:15, cells = 500, balanced = TRUE)

##Determine the 'dimensionality' of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
onefiftyGEX <- JackStraw(onefiftyGEX, num.replicate = 100)
onefiftyGEX <- ScoreJackStraw(onefiftyGEX, dims = 1:20)

JackStrawPlot(onefiftyGEX, dims = 1:15)

ElbowPlot(onefiftyGEX)

##Cluster the cells

onefiftyGEX <- FindNeighbors(onefiftyGEX, dims = 1:10)
onefiftyGEX <- FindClusters(onefiftyGEX, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(onefiftyGEX), 5)

##Run non-linear dimensional reduction (UMAP/tsne)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
onefiftyGEX <- RunUMAP(onefiftyGEX, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(onefiftyGEX, reduction = "umap")

saveRDS(onefiftyGEX, file = "C:/Users/User/Desktop/Seurat Output/")

##Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 1
cluster1.markers <- FindMarkers(onefiftyGEX, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(onefiftyGEX, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
onefiftyGEX.markers <- FindAllMarkers(onefiftyGEX, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
onefiftyGEX.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(onefiftyGEX, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(onefiftyGEX, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(onefiftyGEX, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(onefiftyGEX, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

top10 <- onefiftyGEX.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(onefiftyGEX, features = top10$gene) + NoLegend()

##Assigning cell type identity to clusters

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet", "10")
names(new.cluster.ids) <- levels(onefiftyGEX)
onefiftyGEX <- RenameIdents(onefiftyGEX, new.cluster.ids)
DimPlot(onefiftyGEX, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(onefiftyGEX, file = "C:/Users/User/Desktop/Seurat Output/")
