# Marker Gene Identification for selected inhibitory neuron subtypes
From the snRNA dataset, we divided the inhibitory neuron population into 54 subtypes using [Label Transfer](https://satijalab.org/seurat/archive/v3.0/integration.html) from Seurat. 
This is the list of marker genes for a selected subset of interested subtypes.

For the marker gene identification, we used [FindMarkers()](https://satijalab.org/seurat/reference/findmarkers) from Seurat. You may refer to the script we used to generate the marker gene lists here: `../scripts/5_inh_subtypes/explore_inh_subtypes.ipynb`

Below, we list the associated files in this directory. 

## Marker gene list for selected Inh subtypes
- Significant marker genes (adjusted p-value < 0.05) for each of selected Inh subtypes: `Annotate_Inh_Subtype/Inh_subset_cells/Marker_Genes.xlsx`
