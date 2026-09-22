# Results and interpretation

These findings summarize the supplied project report and its figures. The underlying statistical tables still need to be added and reconciled before a complete release.

## Differential expression

The report records 21,179 retained features and 1,780 significant features for Cd versus Control: 845 upregulated and 935 downregulated at adjusted P < 0.05 and |log₂ fold change| ≥ 1. The two directional groups sum to the reported total. This arithmetic check does not independently validate the statistical analysis.

## Sample relationships

PC1 and PC2 account for 86% and 6% of the variance represented in the PCA input. Cd_2 and Cd_3 are displaced from controls along PC1; Cd_1 is intermediate. In the displayed distance matrix, Cd_1 is closer to CK_3 (37.1) than to Cd_2 (43.3) or Cd_3 (50.0). Its correlation with CK_3 is approximately 0.99, compared with 0.97 with Cd_3. Thus, treatment-associated variation is present, but clustering is not uniform across all replicates.

High overall correlation is compatible with many shared expression patterns and does not by itself verify treatment labels, rule out batch effects, or prove replicate quality. The figures alone do not justify excluding Cd_1. Metadata, alignment/assignment metrics and experimental context should be considered together.

In the top-variable-gene heatmap, colors are relative within each gene. Red indicates a value above that gene's mean across samples and blue a value below it. The color is not an absolute transcript count or the log₂ fold change from DESeq2.

## Functional themes

The report describes enrichment of photosynthesis, photosystem and chloroplast functions among upregulated genes. Other reported themes include redox processes, amino-acid metabolism, phenylpropanoid biosynthesis, lipid and oxylipin pathways, and defense-associated processes. These are hypotheses about transcriptional responses. RNA-seq and enrichment do not establish changes in photosynthetic rate, ROS concentration or enzyme activity.

Fold enrichment compares the fraction of selected genes assigned to a term with the corresponding background fraction. It is not an expression fold change. Full GO/KEGG tables, including adjusted P values and contributing gene IDs, are required to evaluate each term. The GO figure labelled 'recognition of pollen' should not be interpreted as evidence of pollen biology in shoot tissue without inspecting its contributing genes.

## Scope

There are three biological replicates per condition, no new experimental validation, and annotation-dependent enrichment. This is an independent reanalysis of the published experiment, not an independent biological replication. Exact numerical agreement with the original paper is not expected when methods differ. Agreement in broad themes also does not prove agreement for every gene.
