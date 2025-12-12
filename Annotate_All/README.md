# Marker Gene Identification for all cell types in the snRNA dataset
This is the list of marker genes for all cell types (astrocytes, endothelial cells, excitatory neurons, inhibitory neurons, microglia, border-associated macrophage, oligodendrocytes, oligodendrocyte precursor cells and pericytes) using the snRNA dataset.

For the marker gene identification, we used [FindMarkers()](https://satijalab.org/seurat/reference/findmarkers) from Seurat. You may refer to the script we used to generate the marker gene lists here: `scripts/1_all_cell_types/explore_all_ct_rna.ipynb`

Below, we list the associated files in this directory. 

## Marker gene list for all 9 cell types
- Significant marker genes (adjusted p-value < 0.05) for each of the 9 cell types: `significant_celltype_markers.csv`

We also took all cells from the snRNA dataset, performed clustering using Seurat, obtained 36 clusters, and reran marker gene identification for each cluster. 
## Marker gene list for 36 clusters obtained from Seurat clustering
- Significant marker genes (adjusted p-value < 0.05) for each of the 36 clusters from Seurat: `36_clusters/Marker_Genes.xlsx`
