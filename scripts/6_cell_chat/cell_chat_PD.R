library(CellChat)
library(patchwork)
library(Seurat)
library(dplyr)
options(stringsAsFactors = FALSE)

ifnb = readRDS("hc_rbd_pd_joint_integrated.rds")
# already normalized, feature selected, and scaled

ifnb <- subset(ifnb, subset = disease == "PD")

# Rename cell types in the 'celltype' column
ifnb@meta.data$celltype <- recode(ifnb@meta.data$celltype,
  "astrocytes" = "Ast",
  "endothelia cells" = "Endo",
  "excitatory neurons" = "ExN",
  "inhibitory neurons" = "InN",
  "microglia" = "Mic",
  "oligodendrocyte precursor cells" = "OPC",
  "oligodendrocytes" = "Oligo",
  "pericytes" = "Peri"
)

ifnb <- subset(ifnb, subset = celltype != "Peri")

table(ifnb@meta.data$celltype)

data.input <- GetAssayData(ifnb, assay = "RNA", slot = "data") # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
meta <- data.frame(celltype = ifnb@meta.data$celltype, row.names = rownames(ifnb@meta.data))
print(all(rownames(meta) == colnames(data.input)))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
# cellchat <- createCellChat(object = ifnb, group.by = "celltype", assay = "RNA")

CellChatDB <- CellChatDB.human

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)


# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use


ptm = Sys.time()
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 13.20763
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


options(future.globals.maxSize = 2 * 1024^3)  # 2 GiB
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

saveRDS(cellchat, file = "cellchat_popT_PD1.rds")

