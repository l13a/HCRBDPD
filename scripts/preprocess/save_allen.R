# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)

# Step 1: Load Allen Brain Institute Data
allen_matrix <- read.csv("allen_matrix.csv", row.names = 1)
print("read allen matrix")
print(head(rownames(allen_matrix)))
allen_metadata <- read.csv("allen_metadata.csv")
print("read allen metadata")
print(head(allen_metadata))

# Step 2: Preprocess Allen Data
# Ensure rownames in allen_matrix correspond to gene names and columns to cell IDs
allen_matrix <- as.matrix(allen_matrix)
allen_matrix <- t(allen_matrix)

# Add metadata to the Allen data
allen_seurat <- CreateSeuratObject(counts = allen_matrix)
print("sucessfully created seurat object")

rownames(allen_metadata) <- allen_metadata$sample_name
# Subset metadata to include only cells in the Seurat object
allen_metadata <- allen_metadata[colnames(allen_seurat), ]

allen_seurat <- AddMetaData(allen_seurat, metadata = allen_metadata)
print(allen_seurat)

saveRDS(allen_seurat, file = "allen_seurat.rds")