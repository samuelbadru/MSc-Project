# Load required libraries
library(dplyr)
library(Seurat)
library(patchwork)

# Setup the Seurat Object  - Load the dataset and Initialize the Seurat object with the raw (non-normalized data)
your.data <- Read10X(data.dir = "path.to/filtered_feature_bc_matrix/")
your <- CreateSeuratObject(counts = your.data$`Gene Expression`, project = "name.it")

# Create your own clusters: download cvs from loupe browser:
cell_types <- read.csv("path.to/134.csv")

# read what are the column names in your cvs and write in code accordingly (highlighted in yellow): 
clusters <- colnames(your) %>% tibble::enframe(value = "Barcode") %>% left_join(cell_types) %>% pull(`Cell.types`)

#  Add this info to the object using AddMetaData function: https://rdrr.io/github/satijalab/seurat/man/AddMetaData.html
your <- AddMetaData(your, metadata = clusters, col.name = "clusters")
Idents(your) <- 'clusters'

# From here follow up with the "Standard pre-processing workflow" OR

# Alternatively you can add an extra line using the merge function from Seurat to put all your objects into one, thus working only once, rather than repeating this 3 or more times: it's like cellranger aggregate, expect you can merge any kind of your: https://satijalab.org/seurat/v3.1/merge_vignette.html
pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")

# Then you follow up  with "Standard pre-processing workflow"


# 1-50
library(dplyr)
library(Seurat)
library(patchwork)

onefifty.data <- Read10X(data.dir = "C:/Users/User/Desktop/Cellranger Outputs/count/1_50GEXFB/filtered_feature_bc_matrix1_50/Zipped/")
onefiftyGEX <- CreateSeuratObject(counts = onefifty.data[["Gene Expression"]], project = "1-50", min.cells = 3, min.features = 200)
cell_types150 <- read.csv("C:/Users/User/Desktop/Loupe CSV/1-50.csv")
clusters150 <- colnames(onefiftyGEX) %>% tibble::enframe(value = "Barcode") %>% left_join(cell_types150) %>% pull(`Cell.Type`)
onefiftyGEX <- AddMetaData(onefiftyGEX, metadata = clusters150, col.name = "clusters")
Idents(onefiftyGEX) <- 'clusters'

# 2-8
library(dplyr)
library(Seurat)
library(patchwork)

twoeight.data <- Read10X(data.dir = "C:/Users/User/Desktop/Cellranger Outputs/count/2_8GEXFB/filtered_feature_bc_matrix2_8/")
twoeightGEX <- CreateSeuratObject(counts = twoeight.data[["Gene Expression"]], project = "2-8", min.cells = 3, min.features = 200)
cell_types28 <- read.csv("C:/Users/User/Desktop/Loupe CSV/2-8.csv")
clusters28 <- colnames(twoeightGEX) %>% tibble::enframe(value = "Barcode") %>% left_join(cell_types28) %>% pull(`Cell.Type`)
twoeightGEX <- AddMetaData(twoeightGEX, metadata = clusters28, col.name = "clusters")
Idents(twoeightGEX) <- 'clusters'

# 47
library(dplyr)
library(Seurat)
library(patchwork)

fortyseven.data <- Read10X(data.dir = "C:/Users/User/Desktop/Cellranger Outputs/count/47_49GEXFB/filtered_feature_bc_matrix/")
fortysevenGEX <- CreateSeuratObject(counts = fortyseven.data[["Gene Expression"]], project = "47-49", min.cells = 3, min.features = 200)
cell_types47 <- read.csv("C:/Users/User/Desktop/Loupe CSV/47.csv")
clusters47 <- colnames(fortysevenGEX) %>% tibble::enframe(value = "Barcode") %>% left_join(cell_types47) %>% pull(`Cell.Type`)
fortysevenGEX <- AddMetaData(fortysevenGEX, metadata = clusters47, col.name = "clusters")
Idents(fortysevenGEX) <- 'clusters'

# 49
library(dplyr)
library(Seurat)
library(patchwork)

fortynine.data <- Read10X(data.dir = "C:/Users/User/Desktop/Cellranger Outputs/count/49GEX/filtered_feature_bc_matrix/")
fortynine <- CreateSeuratObject(counts = fortynine.data, project = "49", min.cells = 3, min.features = 200)
cell_types49 <- read.csv("C:/Users/User/Desktop/Loupe CSV/49.csv")
clusters49 <- colnames(fortynine) %>% tibble::enframe(value = "Barcode") %>% left_join(cell_types49) %>% pull(`Cell.Type`)
fortynine <- AddMetaData(fortynine, metadata = clusters49, col.name = "clusters")
Idents(fortynine) <- 'clusters'

allsamples <- merge(onefiftyGEX, y = c(twoeightGEX, fortysevenGEX, fortynine), add.cell.ids = c("AZ", "BC", "X", "Z"), project = "Report")