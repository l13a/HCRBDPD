# Human Brain Cell State Dynamics in Prodromal and Early Parkinson’s Disease Progression 
This is the source code repository for *Human Brain Cell State Dynamics in Prodromal and Early Parkinson’s Disease Progression*. The main analyses use R, with specific packages detailed below.

![Analysis Workflow](data/workflow.svg)

## System Requirements

The main language required to replicate the analysis pipeline is R. 
Python is also required to reproduce proportion test results and TF-gene networks. We used:

- R 4.3.0
- Python 3.9.23
  
in our analyses.

Please use the attached `environment.yml` for the full list of system requirements and versions. 

## Installation

To download and install the pipeline, simply clone this repository:

```
git clone https://github.com/l13a/HCRBDPD.git
```

and create a virtual environment with `conda` using the supported `environment.yml`:

```
conda env create -f environment.yml
``` 

## Directory structure

```
.
├── Annotate_All              # Marker genes for all cell types
├── Annotate_Exc_Subtype      # Marker genes for Exc neuron subtypes
├── Annotate_Inh_Subtype      # Marker genes for Inh neurons subtypes
├── Annotate_Mic_Subtype      # Marker genes for Mic subtypes
├── DEG_GO_Analysis           # DEGs and GO results
├── data                      # intermediate saved data 
├── scripts                   # notebooks to reproduce all analysis, divded into sections from the paper
├── environment.yaml          # Reproducible environment via conda
├── LICENSE
└── README.md
```


## Cite Our Work

To be added...
