library(Seurat)
library(Signac)

# Load previously saved objects
hc  <- readRDS("hc_processed_rna+atac.rds")
rbd <- readRDS("rbd_processed_rna+atac.rds")
pd  <- readRDS("pd_processed_rna+atac.rds")

# Helper to keep only ATAC assay and tag disease
prep <- function(obj, label) {
  DefaultAssay(obj) <- "ATAC"
  if ("RNA" %in% Assays(obj)) obj[["RNA"]] <- NULL
  obj$disease <- label
  obj <- DietSeurat(
    obj,
    assays     = "ATAC",
    counts     = TRUE,
    data       = FALSE,
    scale.data = FALSE,
    dimreducs  = NULL,
    graphs     = NULL
  )
  obj
}

hc  <- prep(hc,  "HC")
rbd <- prep(rbd, "RBD")
pd  <- prep(pd,  "PD")

# Merge all three (add.cell.ids is optional; remove if you already have prefixes)
combined <- merge(
  hc,
  y = list(rbd, pd),
  project = "hc_rbd_pd_atac_only"
)
DefaultAssay(combined) <- "ATAC"
combined$disease <- factor(combined$disease, levels = c("HC","RBD","PD"))

# Quick checks
print(Assays(combined))        # should show only "ATAC"
print(table(combined$disease)) # counts per disease

# Save ATAC samples from all disease status to rds
saveRDS(combined, "hc_rbd_pd_processed_atac.rds")