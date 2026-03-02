#!/usr/bin/env Rscript


###############
## Input files ##
###############

#
#
#
#

###############
## Run Pipeline ##
###############

## NB : replace files and folders paths before running (top documents)

library(rmarkdown)

### Pre-processing

# Pre-processing of omics and DSRT data (R v3.6.0, MOFA2 v0.99.7)

## Merge cells markers (R v4.4.2)
render("~/Code/MOFA_AML/src/3_Analysis_local/1_processing/pre_processing_AML_cell_types_markers.Rmd")

### Run MOFA (R v3.6.0, MOFA2 v0.99.7)
render("~/Code/MOFA_AML/src/1_MOFA_training/MOFA_v04_train.r")
render("~/Code/MOFA_AML/src/2_Analysis_server/Downstream_analysis_MOFA.Rmd")

### Validation factors

# CITE-seq data based scores 
render("~/Code/MOFA_AML/src/2_Analysis_server/Clinseq_pro_RNAseq_AML_cell_types_vs_factors.Rmd") #(R v3.6.0, Seurat v3.2.2)
render("~/Code/MOFA_AML/src/3_Analysis_local/2_validation/CITE_seq_association_scores_cell_types_surface_markers.Rmd") #(R v4.4.2, Seurat v5.3)

# Transcriptomics based scores and survival analysis (R v4.4.2)
render("~/Code/MOFA_AML/src/3_Analysis_local/2_validation/hallmarks_scores_based_transcriptomics_survival.Rmd")

### Figures (R v4.4.2)
render("~/Code/MOFA_AML/src/3_Analysis_local/3_figures/extract_top_drugs_rank_factors_MOFA.Rmd")
render("~/Code/MOFA_AML/src/3_Analysis_local/3_figures/figures_MOFA_variance_per_factor.Rmd")
render("~/Code/MOFA_AML/src/3_Analysis_local/3_figures/heatmap_half_triangles_pathways_MOFA.Rmd")
render("~/Code/MOFA_AML/src/3_Analysis_local/3_figures/sankey_plot_cell_types_MOFA.Rmd")

