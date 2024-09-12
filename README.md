TITLE PAPER
==========================

This repository contains analysis code for the RNA-seq expression project on AD mouse model at two different time points. This project was carried out by researchers at the [Wang Lab, MUSC](https://medicine.musc.edu/departments/neuroscience/research/cowan) and [Berto Lab, MUSC](https://bertolab.org/).

## Cite this

If you use anything in this repository please cite the following publication:

URL: 


## Files

| directory | contents | code |
| --------- | -------- | -------- |
| [`futcounts`](futcounts/) | Input data of the initial processing and quality check. | 01_Create_InputData.R|
| [`dge`](dge/) | Output of the Differential expression analysis. | 02_Dge.R \ 03_Convert_ID.R \ 04_DGE_Viz.R \ 05_DGE_GeneOntology.R|
| [`utils`](utils/) | Utility functions and data. | 08_Enrichment_Modules.sh |
