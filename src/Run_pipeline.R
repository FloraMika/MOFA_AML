#!/usr/bin/env Rscript

###############
## Input files ##
###############

load("/home/flomik/Data/MOFA_R_Objects/AML_MAexp_v05.RData")
load("/home/flomik/Data/MOFA_R_Objects/MOFA_v04_train_list.RData")

# AML_MAexp_v05.RData : MultiAssayExperiment object with filtered omics and DSRT data
# MOFA_v04_train_list.RData : MOFA output including MOFA model
# Wells_combined_sctransformed.rds : Seurat object with CITE-seq data
# Ensembl2Reactome_All_Levels.txt : GMT file including reactome pathways and associated genes
# Clinical_mofa_retro.xlsx : dataframe including clinal parameters for retrospective cohort

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

### Review (R v4.4.2)

#### extra figures
render("~/Code/MOFA_AML/src/4_review/1_figures/Tojo/hallmark_circularplot.r")
render("~/Code/MOFA_AML/src/4_review/1_figures/1_review_figures.Rmd")
#### Robustness
render("~/Code/MOFA_AML/src/4_review/2_robustness/1_MOFA_robustness_generate_data_run_models.Rmd")
render("~/Code/MOFA_AML/src/4_review/2_robustness/2_MOFA_robustness_generate_data_run_models_part_2.Rmd")
render("~/Code/MOFA_AML/src/4_review/2_robustness/3_MOFA_robustness_correlations_2.Rmd")

#### validation in single cell
render("~/Code/MOFA_AML/src/4_review/3_scRNA/1_KI_scRNA_analysis_BMMap_v2.rmd")
render("~/Code/MOFA_AML/src/4_review/3_scRNA/2_Validation_scRNA_analysis_BMMap_Zheng_v2.rmd")

#### change sign
render("~/Code/MOFA_AML/src/4_review/4_sign/1_change_sign_figures_single_cells.rmd")
render( "~/Code/MOFA_AML/src/4_review/4_sign/1_change_sign_figures_single_cells.rmd")