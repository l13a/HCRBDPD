library(Signac)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

##### Combine raw RNA + ATAC data files + metadata into rds object for PD samples

df <- read.csv("HCTXJ_CTR_PFC.csv", stringsAsFactors = FALSE)
df <- df %>% dplyr::filter(library_type == "Gene Expression")
row_idxes <- which(grepl("^PD", df$sample))
# keep track of all batch names
sample_names <- df$sample[row_idxes]
ids <- row_idxes - 1
print(row_idxes)
print(sample_names)

id = ids[1]
print(paste0("For object number ", 1, " , with id ", id, " , and row idx ", row_idxes[1], " , and batch ", df$sample[row_idxes[1]]))
rna = LoadH5Seurat(paste0("processed_out/test_rna_id", id, ".h5Seurat"))
celltype = read.csv(paste0("ct_list/celltype_id", id, ".csv"))    
rna$celltype = celltype$ct
data_list = list(rna)
rna$batch = df$sample[row_idxes[1]]

# only need first 9 objects for HC-RNA seurat since contains ATAC assay
for (i in 2:9){
    id = ids[i]
    print(paste0("For object number ", i, " , with id ", id, " , and row idx ", row_idxes[i], " , and batch ", df$sample[row_idxes[i]]))
    rna = LoadH5Seurat(paste0("processed_out/test_rna_id", id, ".h5Seurat"))
    celltype = read.csv(paste0("ct_list/celltype_id", id, ".csv"))    
    rna$celltype = celltype$ct
    rna$batch = df$sample[row_idxes[i]] # store corresponding batch info for id i
    data_list = append(data_list,rna)    
}

print(data_list[[1]])
print(data_list[[2]])

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

for (i in 1:9){
    id = ids[i]
    fragpath = paste0("paireddata_id", id, "/outs/atac_fragments.tsv.gz")
    print(fragpath)
    # create ATAC assay and add it to the object
    data_list[[i]][["ATAC"]] <- CreateChromatinAssay(
      counts = data_list[[i]][["ATAC"]]$counts,
      sep = c(":", "-"),
      fragments = fragpath,
      annotation = annotation
    )
}

ifnb <-  merge(data_list[[1]], y = c(data_list[[2]], 
                                     data_list[[3]], 
                                     data_list[[4]], 
                                     data_list[[5]], 
                                     data_list[[6]], 
                                     data_list[[7]], 
                                     data_list[[8]], 
                                     data_list[[9]]
                                    ), add.cell.ids = sample_names, project = "pd_rna_alldata")

saveRDS(ifnb, file = "pd_rna+atac.rds")
print("Done saving!!")