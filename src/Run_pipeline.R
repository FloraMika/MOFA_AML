#!/usr/bin/env Rscript

###############
## Input files ##
###############

# AML_MAexp_v05.RData : MultiAssayExperiment object with omics and DSRT data
# pData : dataframe including clinal parameters for retrospective cohort
# AML_MAexp_v05_filt.RData : MultiAssayExperiment object with filtered and scaled omics and DSRT data
# MOFA_v04_train_list.RData : MOFA output including MOFA model
# Wells_combined_sctransformed.rds : Seurat object with CITE-seq data
# Ensembl2Reactome_All_Levels.txt : GMT file including reactome pathways and associated genes


###############
## Run Pipeline ##
###############

## put data input folder into MOFA_AML/data

### set Path
path_general <- getwd()
new_path_general <- gsub("MOFA_AML/.*","MOFA_AML/src/",path)

library(rmarkdown)

### Pre-processing

# Pre-processing of omics and DSRT data (R v3.6.0, MOFA2 v0.99.7)

### Run MOFA (R v3.6.0, MOFA2 v0.99.7)
render(paste0(new_path_general, "1_MOFA_training/1_MOFA_v04_train.Rmd"))
render(paste0(new_path_general, "1_MOFA_training/2_Downstream_analysis_MOFA.Rmd"))

# Transcriptomics based scores and survival analysis (R v4.4.2)
render(paste0(new_path_general, "2_Analysis_local/1_survival/hallmarks_scores_based_transcriptomics_survival_v2.Rmd"))

### Figures (R v4.4.2)
render(paste0(new_path_general, "2_Analysis_local/2_figures/1_figures_MOFA_variance_per_factor.Rmd"))
render(paste0(new_path_general, "2_Analysis_local/2_figures/2_heatmap_half_triangles_pathways_MOFA.Rmd"))
render(paste0(new_path_general, "2_Analysis_local/2_figures/3_extract_top_drugs_rank_factors_MOFA.Rmd"))
render(paste0(new_path_general, "2_Analysis_local/2_figures/4_scatter_plots_drugs_factors.Rmd"))
render(paste0(new_path_general, "2_Analysis_local/2_figures/hallmark_circularplot_Tojo.r"))


### Review (R v4.4.2)

#### extra figures
render(paste0(new_path_general, "3_review/1_figures/1_review_figures.Rmd"))

#### Robustness
render(paste0(new_path_general, "3_review/2_robustness/1_MOFA_robustness_generate_data_run_models.Rmd"))
render(paste0(new_path_general, "3_review/2_robustness/2_MOFA_robustness_generate_data_run_models_part_2.Rmd"))
render(paste0(new_path_general, "3_review/2_robustness/3_MOFA_robustness_correlations_2.Rmd"))

#### validation in single cell
render(paste0(new_path_general, "3_review/3_scRNA/1_KI_scRNA_analysis_BMMap_v2.rmd"))
render(paste0(new_path_general, "3_review/3_scRNA/2_Validation_scRNA_analysis_BMMap_Zheng_v2.rmd"))

#### change sign
render(paste0(new_path_general, "3_review/4_sign/2_change_sign_figures_beatAML.rmd"))

