# MOFA_AML

Code for the analyses from the paper : Data-driven hallmarks in acute myeloid leukemia: a strategy towards diagnostics and individualized therapy 

## About

In this study, we introduce the Data-Driven Hallmark (DDHM), a low dimensional representation concept for facilitating multi-omics interpretation in translational cancer research, using acute myeloid leukemia (AML) as a model. MOFA was used for integration of multi-omics data (genomics, methylation, transcriptomics, proteomics) and ex vivo drug response profiles from 118 AML patients at diagnosis, revealing complex disease related mechanisms of AML.

## Installation

## Clone repository
```
git clone https://github.com/FloraMika/MOFA_AML.git
cd MOFA_AML
```

## Requirements

1. A linux distribution.
2. Install the few required **R packages** :

```
Rscript source/requierements.R
# This command will install the following packages:
# sva==3.26.0
# pca3d==0.10
# dplyr==0.8.0.1
```

## Pipeline

![Pipeline](img/#14_venn_diagramm.png)


## AML MOFA model

MOFA model can be accessed using the following code. For the installation and use of MOFA2, we recommend looking at the authors website https://biofam.github.io/MOFA2/tutorials.html.

Analyses can be reproduced using the different R codes available in the src folder.

![Top weighted features for each data type for 11 hallmarks](img/#14_venn_diagramm.png)

## Authors
- [Flora Mikaeloff](https://github.com/FloraMika)
- Tojo James
- Francesco Marabita

## Acknowledgment

## License