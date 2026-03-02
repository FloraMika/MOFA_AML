#!/usr/bin/env Rscript

###############
## Run Pipeline ##
###############

library(rmarkdown)

### Training MOFA model (R v3.6.0, MOFA2 v0.99.7)

# Pre-processing of omics and DSRT data
# TO ADD
# Model training
render("src/1_MOFA_training/MOFA_v04_train.r")


### Downstream analysis (R v3.6.0, MOFA2 v0.99.7)

# Downstream analysis
render("src/2_Analysis_server/Downstream_analysis_MOFA.Rmd")


### Validation factors

# Merge markers cell types (R v4.4.2)
render("src/3_Analysis_local/1_processing/pre_processing_AML_cell_types_markers.Rmd")

# CITE-seq data based scores (R v4.4.2, Seurat v3.2.2)
render("src/2_Analysis_server/Clinseq_pro_RNAseq_AML_cell_types_vs_factors.Rmd")
render("src/3_Analysis_local/2_validation/CITE_seq_association_scores_cell_types_surface_markers.Rmd")

# Transcriptomics based scores and survival analysis (R v4.4.2)
render("src/3_Analysis_local/2_validation/hallmarks_scores_based_transcriptomics_survival.Rmd")

### Figures (R v4.4.2)

render("src/3_Analysis_local/3_figures/extract_top_drugs_rank_factors_MOFA.Rmd")
render("src/3_Analysis_local/3_figures/figures_MOFA_variance_per_factor.Rmd")
render("src/3_Analysis_local/3_figures/heatmap_half_triangles_pathways_MOFA.Rmd")
render("src/3_Analysis_local/3_figures/sankey_plot_cell_types_MOFA.Rmd")

