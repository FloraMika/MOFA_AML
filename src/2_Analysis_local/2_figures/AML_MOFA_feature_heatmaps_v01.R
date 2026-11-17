library(MultiAssayExperiment)
library(SummarizedExperiment)
library(RColorBrewer)
library(MOFA2)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(patchwork)
library(minfi)
library(tidyverse)

load("~/Projects/AML/out/MOFA_v04/AML_MOFAplus_sel_model_v04.RData")
load("~/Projects/AML/data/DNAm/AML_DNAm_GenomicRatioSet_v01.RData")

# Factor 1
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F1_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor1", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor1", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor1", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor1=value) %>% 
  arrange(Factor1) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor1", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor1", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor1", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor1=value) %>% 
  arrange(Factor1) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor1", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor1", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor1", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor1=value) %>% 
  arrange(Factor1) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.DNAm

## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor1", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor1", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor1=value) %>%
  arrange(Factor1) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order =  features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint

## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor1= circlize::colorRamp2(c(min(annot.table$Factor1), 0, max(annot.table$Factor1)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()



# Factor 2
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F2_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor2", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor2", nfeatures = 20, sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor2", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor2=value) %>% 
  arrange(Factor2) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
        ) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor2", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor2", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor2", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor2=value) %>% 
  arrange(Factor2) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
        ) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor2", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor2", nfeatures = 20, sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor2", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor2=value) %>% 
  arrange(Factor2) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
        ) -> ht.DNAm


## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor2", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor2", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor2=value) %>%
  arrange(Factor2) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order = features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint

## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor2= circlize::colorRamp2(c(min(annot.table$Factor2), 0, max(annot.table$Factor2)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()

pdf("~/Projects/AML/figures/AML_MOFA_feature_weights_F2_v01.pdf", width = 7, height = 7)
options(ggrepel.max.overlaps = Inf)
plot_top_weights(mofa.sel,view = "DSRT", factor = "Factor2", nfeatures = 10, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DSRT", factor = "Factor2", nfeatures = 10, sign = "negative") -> p1.neg

plot_weights(mofa.sel, view = "DSRT", factors ="Factor2",
             manual = c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),
             color_manual = colorRampPalette(c("red", "blue"))(20)) -> p1.weights
p1.weights$layers[[2]]$data$feature <- features_metadata(mofa.sel)$name[match(p1.weights$layers[[2]]$data$feature,features_metadata(mofa.sel)$feature)]
p1.weights +
  scale_y_discrete(expand = expansion(mult = 0.1))
dev.off()

# Factor 3
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F3_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor3", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor3", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>% 
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor3", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor3", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>% 
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor3", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor3", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>% 
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.DNAm

## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor3", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>%
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order =  features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint

## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor3= circlize::colorRamp2(c(min(annot.table$Factor3), 0, max(annot.table$Factor3)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()

# Factor 5
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F5_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor5", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor5", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor5", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor5=value) %>% 
  arrange(Factor5) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor5", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor5", nfeatures = 20, sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor5", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor5=value) %>% 
  arrange(Factor5) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor5", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor5", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor5", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor5=value) %>% 
  arrange(Factor5) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.DNAm

## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor5", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor5=value) %>%
  arrange(Factor5) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order =  features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint

## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor5= circlize::colorRamp2(c(min(annot.table$Factor5), 0, max(annot.table$Factor5)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()


# Factor 8
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F3_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor3", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor3", nfeatures = 20, sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>% 
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor3", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor3", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>% 
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor3", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor3", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>% 
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.DNAm

## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor3", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor3", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor3=value) %>%
  arrange(Factor3) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order = features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint


## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor3= circlize::colorRamp2(c(min(annot.table$Factor3), 0, max(annot.table$Factor3)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()

# Factor 8
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F8_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor8", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor8", nfeatures = 20, sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor8", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor8=value) %>% 
  arrange(Factor8) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor8", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor8", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor8", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor8=value) %>% 
  arrange(Factor8) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor8", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor8", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor8", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor8=value) %>% 
  arrange(Factor8) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.DNAm

## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor8", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor8", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor8=value) %>%
  arrange(Factor8) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order = features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint


## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor8= circlize::colorRamp2(c(min(annot.table$Factor8), 0, max(annot.table$Factor8)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()

# Factor 9
pdf("~/Projects/AML/figures/AML_MOFA_feature_heatmaps_F9_v01.pdf", width = 10, height = 13)

## Transcripts
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor9", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "RNAseq", factor = "Factor9", nfeatures = 20, sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "RNAseq")$RNAseq$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor9", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor9=value) %>% 
  arrange(Factor9) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Transcript Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.rna

## Proteins
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor9", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "Proteomics", factor = "Factor9", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "Proteomics")$Proteomics$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- mapIds(org.Hs.eg.db, keys = sapply(strsplit(rownames(toPlot), ".", fixed = T),"[", 1), column = "SYMBOL", keytype = "ENSEMBL")

get_factors(mofa.sel,factors = "Factor9", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor9=value) %>% 
  arrange(Factor9) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(t(scale(t(toPlot))), name = "Protein Z-score",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.prot

## DNAm
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor9", nfeatures = 20, sign = "positive") -> p1.pos
plot_top_weights(mofa.sel,view = "DNAm", factor = "Factor9", nfeatures = 20,sign = "negative") -> p1.neg

toPlot <- get_data(mofa.sel, views = "DNAm")$DNAm$AML
toPlot <- toPlot[c(as.character(p1.pos$data$feature), rev(as.character(p1.neg$data$feature))),]
rownames(toPlot) <- paste0(as.character(seqnames(data.Rset[rownames(toPlot),])),
                           ":",
                           end(data.Rset[rownames(toPlot),]),
                           " ",
                           features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)])

get_factors(mofa.sel,factors = "Factor9", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor9=value) %>% 
  arrange(Factor9) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>% 
  column_to_rownames("ID") -> annot.table

Heatmap(toPlot, name = "DNAm M-value",
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F,
        row_names_gp = gpar(fontsize = 8)
) -> ht.DNAm

## Mutations oncoprint
plot_top_weights(mofa.sel,view = "Mutations", factor = "Factor9", nfeatures = Inf, sign="all") -> p1.mut

toPlot <- get_data(mofa.sel, views = "Mutations")$Mutations$AML
rownames(toPlot) <- features_metadata(mofa.sel)$name[match(rownames(toPlot), features_metadata(mofa.sel)$id)]

get_factors(mofa.sel,factors = "Factor9", as.data.frame = T) %>%
  select(sample, value) %>%
  rename(ID=sample, Factor9=value) %>%
  arrange(Factor9) -> annot.table

toPlot <- toPlot[,annot.table$ID]

annot.table[match(colnames(toPlot), annot.table$ID),] %>%
  remove_rownames() %>%
  column_to_rownames("ID") -> annot.table

tmp <- toPlot
tmp[tmp==1] <- "Presence"
tmp[tmp==0] <- ""

col=c(Presence="blue")

oncoPrint(tmp,
          get_type = function(x) {x},
          alter_fun = list(
            Presence = function(x, y, w, h) grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = col["Presence"], col = NA))
          ),
          col=col, show_pct = F, top_annotation = NULL, heatmap_legend_param = list(title="Mutations"),
          row_names_gp = gpar(fontsize = 8),
          row_order = features_metadata(mofa.sel)$name[match(as.character(p1.mut$data$feature),features_metadata(mofa.sel)$id)]) -> ht.oncoprint


## top annotation
top_annotation = HeatmapAnnotation(df=annot.table,
                                   col = list(Factor9= circlize::colorRamp2(c(min(annot.table$Factor9), 0, max(annot.table$Factor9)), c("blue", "yellow", "red"))))

# draw heatmaps
top_annotation %v% ht.rna %v% ht.prot %v% ht.DNAm %v% ht.oncoprint

dev.off()

