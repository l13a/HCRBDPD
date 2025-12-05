library(Seurat)
library(ggplot2)
library(enrichR)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(scales)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(tibble)
library(grid)
library(Matrix)

## Integrate allen brain subtypes to our RNA dataset using label transfer

allen_seurat = readRDS("allen_seurat.rds")
allen_seurat <- NormalizeData(allen_seurat)
allen_seurat <- FindVariableFeatures(allen_seurat)
allen_seurat <- ScaleData(allen_seurat)

# Subset only excitatory neurons
cells_to_keep <- rownames(allen_seurat@meta.data)[grepl("^Exc", allen_seurat@meta.data$cluster_label)]
allen_seurat_Ex <- subset(allen_seurat, cells = cells_to_keep)
allen_seurat_Ex@meta.data$source <- "allen"

rm(allen_seurat)


ifnb = readRDS("hc_rbd_pd_joint_integrated.rds") # already normalized, feature selected, and scaled
ifnb_testct <- subset(ifnb, subset = celltype == "excitatory neurons")
ifnb_testct@meta.data$source <- "Zhang"

rm(ifnb)

ifnb_testct <- FindNeighbors(ifnb_testct, reduction = "integrated.cca", dims = 1:30)
ifnb_testct <- FindClusters(ifnb_testct, resolution = 0.5)
ifnb_testct <- RunUMAP(ifnb_testct, dims = 1:30, reduction = "integrated.cca")

# Integrate
anchors <- FindIntegrationAnchors(object.list = list(ifnb_testct, allen_seurat_Ex), 
                                   anchor.features = intersect(rownames(ifnb_testct), rownames(allen_seurat_Ex)))
print("Done finding anchors")
integrated_data <- IntegrateData(anchorset = anchors)
print(integrated_data)

rm(anchors)

# Perform PCA and clustering on the integrated data
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data)
integrated_data <- FindNeighbors(integrated_data, dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

saveRDS(ifnb_testct, file = "integrated_data.rds")

# Transfer labels from Allen data to your data
transfer_anchors <- FindTransferAnchors(reference = allen_seurat_Ex, query = ifnb_testct, 
                                        dims = 1:30)
print("Done finding transfer anchors")

predicted_labels <- TransferData(anchorset = transfer_anchors, refdata = allen_seurat_Ex@meta.data$cluster_label)
ifnb_testct <- AddMetaData(ifnb_testct, metadata = predicted_labels)

saveRDS(ifnb_testct, file = "integrated_ifnb_testct_Ex.rds")