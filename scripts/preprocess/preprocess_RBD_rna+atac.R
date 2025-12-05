library(Signac)
library(Seurat)
library(SeuratDisk)
library(dplyr)

# RNA + ATAC rds object for RBD
ifnb <- readRDS("rbd_rna+atac.rds")
metadata <- ifnb@meta.data
metadata$batch <- sub("(_MAH).*", "\\1", rownames(metadata))
ifnb@meta.data <- metadata

# filter out low quality cells
ifnb <- subset(
  x = ifnb,
  subset = nCount_ATAC < 10000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 0 &
    nucleosome_signal < 1.8 &
    nucleosome_signal > 0.5 &
    TSS.enrichment > 1 &
    TSS.enrichment < 10
)

options(future.globals.maxSize = 8 * 1024^3)
DefaultAssay(ifnb) <- "RNA"
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)

DefaultAssay(ifnb) <- "ATAC"
ifnb <- FindTopFeatures(ifnb, min.cutoff = 5)
ifnb <- RunTFIDF(ifnb)
ifnb <- RunSVD(ifnb)

saveRDS(ifnb, file = "rbd_processed_rna+atac.rds")
print("Done saving!!")