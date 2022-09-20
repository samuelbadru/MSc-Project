##Setup the Seurat Object

install.packages("pacman")

pacman::p_load(dplyr, Seurat, patchwork)

# Load the PBMC dataset
fortynine.data <- Read10X(data.dir = "C:/Users/User/Desktop/Cellranger Outputs/count/49GEX/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
fortynine <- CreateSeuratObject(counts = fortynine.data, project = "49", min.cells = 3, min.features = 200)
fortynine

##Standard pre-processing workflow

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
fortynine[["percent.mt"]] <- PercentageFeatureSet(fortynine, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(fortynine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(fortynine, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fortynine, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

fortynine <- subset(fortynine, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##Normalising the data

fortynine <- NormalizeData(fortynine, normalization.method = "LogNormalize", scale.factor = 10000)

##Identification of highly variable features (feature selection)

fortynine <- FindVariableFeatures(fortynine, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(fortynine), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(fortynine)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

##Scaling the data

all.genes <- rownames(fortynine)
fortynine <- ScaleData(fortynine, features = all.genes)

##Perform linear dimensional reduction

fortynine <- RunPCA(fortynine, features = VariableFeatures(object = fortynine))

# Examine and visualize PCA results a few different ways
print(fortynine[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(fortynine, dims = 1:2, reduction = "pca")

DimPlot(fortynine, reduction = "pca")

DimHeatmap(fortynine, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(fortynine, dims = 1:15, cells = 500, balanced = TRUE)

##Determine the 'dimensionality' of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
fortynine <- JackStraw(fortynine, num.replicate = 100)
fortynine <- ScoreJackStraw(fortynine, dims = 1:20)

JackStrawPlot(fortynine, dims = 1:15)

ElbowPlot(fortynine)

# Look at cluster IDs of the first 5 cells
head(Idents(fortynine), 5)

##Run non-linear dimensional reduction (UMAP/tsne)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
fortynine <- RunUMAP(fortynine, dims = 1:12)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(fortynine, reduction = "umap", label = TRUE)

saveRDS(fortynine, file = "C:/Users/User/Desktop/Seurat Output/")

##Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 1
clusterMB.markers <- FindMarkers(fortynine, ident.1 = 'MB', min.pct = 0.25)
head(clusterMB.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(fortynine, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
fortynine.markers <- FindAllMarkers(fortynine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fortynine.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(fortynine, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(fortynine, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(fortynine, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(fortynine, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                      "CD8A"))

top10 <- fortynine.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(fortynine, features = top10$gene) + NoLegend()

##Assigning cell type identity to clusters

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet", "10")
names(new.cluster.ids) <- levels(fortynine)
fortynine <- RenameIdents(fortynine, new.cluster.ids)
DimPlot(fortynine, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(fortynine, file = "C:/Users/User/Desktop/Seurat Output/")
