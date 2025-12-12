# Marker Gene Identification for all excitatory neuron RORB subtypes
From the snRNA dataset, we divided the excitatory neuron population into 50 subtypes using [Label Transfer](https://satijalab.org/seurat/archive/v3.0/integration.html) from Seurat. 
This is the list of marker genes for all subtypes relating to RORB cells.

For the marker gene identification, we used [FindMarkers()](https://satijalab.org/seurat/reference/findmarkers) from Seurat. You may refer to the script we used to generate the marker gene lists here: `../scripts/4_exc_subtypes/explore_exc_subtypes.ipynb`

Below, we list the associated files in this directory. 

## Marker gene list for all Exc RORB subtypes
- Significant marker genes (adjusted p-value < 0.05) for each of the Exc RORB subtypes: `Annotate_Exc_Subtype/RORB_cells/Marker_Genes.xlsx`
