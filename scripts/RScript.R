library(DESeq2)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)

# Project folders
project_dir = "D:/Nancy/Project"
count_dir = file.path(project_dir, "05_featurecounts")
result_dir = file.path(project_dir, "06_DESeq2")

plot_dir = file.path(result_dir, "plots")
table_dir = file.path(result_dir, "tables")
normalized_dir = file.path(result_dir, "normalized_counts")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(normalized_dir, showWarnings = FALSE, recursive = TRUE)


# Sample information

metadata_file <- file.path(project_dir, "metadata.csv")

sample_info <- read.csv(
  metadata_file,
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "sample_id",
  "run_accession",
  "FASTQ R1",
  "FASTQ R2",
  "treatment",
  "Cd_concentration_uM",
  "duration_days",
  "tissue",
  "organism",
  "replicate"
)

missing_columns <- setdiff(required_columns, names(sample_info))

if (length(missing_columns) > 0) {
  stop(
    "The following metadata columns are missing: ",
    paste(missing_columns, collapse = ", "),
    "\n\nColumns actually found:\n",
    paste(names(sample_info), collapse = ", ")
  )
}

sample_ids <- trimws(sample_info$sample_id)

if (anyNA(sample_ids) || any(sample_ids == "")) {
  stop("Sample IDs are missing.")
}

if (anyDuplicated(sample_ids)) {
  stop("Sample IDs are duplicated.")
}

if (anyDuplicated(sample_info$run_accession)) {
  stop("Run accessions are duplicated.")
}

treatment_clean <- tolower(trimws(sample_info$treatment))

if (!all(treatment_clean %in% c("control", "cd stress"))) {
  stop("Treatment must be either 'control' or 'Cd stress'.")
}

sample_info$condition <- factor(
  ifelse(treatment_clean == "control", "Control", "Cd"),
  levels = c("Control", "Cd")
)

rownames(sample_info) <- sample_ids

# Save sample metadata

write.csv(
  sample_info,
  file.path(result_dir, "sample_metadata.csv"),
  row.names = FALSE
)


# 2. Read counts

count_file = file.path(count_dir, "gene_counts.txt")

featurecounts_raw = read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

count_columns = c(
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

counts = featurecounts_raw[
  ,
  setdiff(colnames(featurecounts_raw), count_columns),
  drop = FALSE
]

rownames(counts) = featurecounts_raw$Geneid

# Convert count to numeric
counts <- as.data.frame(
  lapply(counts, as.numeric),
  row.names = rownames(counts),
  check.names = FALSE
)

# Remove BAM from featureCounts column names.
count_ids <- sub(
  "\\.sorted\\.bam$|\\.bam$",
  "",
  basename(colnames(counts)),
  ignore.case = TRUE
)

# Convert run accessions, such as SRR9113521,
# to analysis sample IDs, such as CK_1.
run_to_sample <- setNames(
  sample_info$sample_id,
  sample_info$run_accession
)

if (all(count_ids %in% names(run_to_sample))) {
  
  colnames(counts) <- unname(run_to_sample[count_ids])
  
} else if (all(count_ids %in% sample_info$sample_id)) {
  
  # Counts are already named CK_1, CK_2, etc.
  colnames(counts) <- count_ids
  
} else {
  
  stop(
    "Count-matrix names do not match either run_accession or sample_id.\n",
    "Counts: ", paste(count_ids, collapse = ", "), "\n",
    "Run accessions: ",
    paste(sample_info$run_accession, collapse = ", "), "\n",
    "Sample IDs: ",
    paste(sample_info$sample_id, collapse = ", ")
  )
}

if (anyDuplicated(colnames(counts))) {
  stop("Duplicate sample names were produced in the count matrix.")
}

if (!setequal(colnames(counts), sample_ids)) {
  stop("Counts and metadata contain different samples.")
}

# Arrange counts in exactly the same order as the metadata.
counts <- counts[, sample_ids, drop = FALSE]

if (!identical(colnames(counts), rownames(sample_info))) {
  stop("Counts and metadata are not in the same order.")
}

#convert the count matrix to numeric matrix
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

# Filter low-count genes
# Keep genes with >=10 counts in at least 3 samples.
keep = rowSums(counts >= 10) >= 3
counts_filtered = counts[keep, ]

cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering :", nrow(counts_filtered), "\n")

# Save filtered raw counts

write.csv(
  counts_filtered,
  file.path(table_dir, "filtered_raw_counts.csv")
)


# Differential expression with DESeq2

coldata <- data.frame(
  condition = sample_info$condition,
  row.names = sample_ids
)

if (!identical(colnames(counts_filtered), rownames(coldata))) {
  stop("Filtered counts and DESeq2 metadata are not in the same order.")
}

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData = coldata,
  design = ~ condition
)

dds = DESeq(dds)

# Save DESeq2 object

saveRDS(
  dds,
  file.path(result_dir, "dds.rds")
)


# Save DESeq2 size factors

size_factors = data.frame(
  sample_id = names(sizeFactors(dds)),
  size_factor = as.numeric(sizeFactors(dds))
)

write.csv(
  size_factors,
  file.path(table_dir, "DESeq2_size_factors.csv"),
  row.names = FALSE
)


# Cd vs Control
res = results(
  dds,
  contrast = c("condition", "Cd", "Control")
)

res = res[order(res$padj), ]

write.csv(
  as.data.frame(res),
  file.path(table_dir, "DESeq2_Cd_vs_Control_ALL.csv")
)


# Define significant DEGs

# Significance threshold:
# adjusted P < 0.05 and |log2FC| >= 1

sig = res[
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) >= 1,
]

deg_up = sig[sig$log2FoldChange >= 1, ]
deg_down = sig[sig$log2FoldChange <= -1, ]

cat("Significant DEGs:", nrow(sig), "\n")
cat("Upregulated DEGs:", nrow(deg_up), "\n")
cat("Downregulated DEGs:", nrow(deg_down), "\n")

write.csv(
  as.data.frame(sig),
  file.path(table_dir, "DEGs_Cd_vs_Control.csv")
)

write.csv(
  as.data.frame(deg_up),
  file.path(table_dir, "DEGs_UP_Cd_vs_Control.csv")
)

write.csv(
  as.data.frame(deg_down),
  file.path(table_dir, "DEGs_DOWN_Cd_vs_Control.csv")
)


# Normalized counts

normalized_counts = counts(dds, normalized = TRUE)

write.csv(
  normalized_counts,
  file.path(
    normalized_dir,
    "normalized_counts.csv"
  )
)


# MA plot

png(
  file.path(
    plot_dir,
    "MA_plot_Cd_vs_Control.png"
  ),
  width = 1200,
  height = 900,
  res = 150
)

plotMA(
  res,
  ylim = c(-8, 8),
  main = "MA Plot - Cd vs Control"
)

dev.off()


# PCA

vsd = vst(
  dds,
  blind = FALSE
)

# Save VST object

saveRDS(
  vsd,
  file.path(result_dir, "vsd.rds")
)

pca_data = plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance = round(
  100 * attr(pca_data, "percentVar")
)

pca_plot = ggplot(
  pca_data,
  aes(PC1, PC2, label = name)
) +
  geom_point(
    aes(shape = group),
    size = 4
  ) +
  geom_text(
    vjust = -1,
    size = 3.5
  ) +
  xlab(
    paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    )
  ) +
  ggtitle(
    "PCA - Rice Shoot RNA-seq"
  ) +
  theme_classic()

ggsave(
  file.path(
    plot_dir,
    "PCA_Cd_vs_Control.png"
  ),
  pca_plot,
  width = 7,
  height = 6,
  dpi = 300
)


# Sample correlation heatmap

sample_cor = cor(
  assay(vsd)
)

png(
  file.path(
    plot_dir,
    "sample_correlation_heatmap.png"
  ),
  width = 1000,
  height = 1000,
  res = 150
)

pheatmap(
  sample_cor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Sample Correlation"
)

dev.off()

write.csv(
  sample_cor,
  file.path(
    table_dir,
    "sample_correlation.csv"
  )
)


# Sample distance heatmap

sample_distance = dist(
  t(assay(vsd))
)

png(
  file.path(
    plot_dir,
    "sample_distance_heatmap.png"
  ),
  width = 1000,
  height = 1000,
  res = 150
)

pheatmap(
  as.matrix(sample_distance),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Sample Distance"
)

dev.off()

write.csv(
  as.matrix(sample_distance),
  file.path(
    table_dir,
    "sample_distance_matrix.csv"
  )
)


# Top 50 variable genes

gene_variance = apply(
  assay(vsd),
  1,
  var
)

top50 = names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[1:min(50, length(gene_variance))]

#This creates a smaller matrix containing only the top 50 variable genes.

mat_top50 = assay(vsd)[
  top50,
  ,
  drop = FALSE
]

# Scale each gene across samples

mat_top50_scaled = t(
  scale(
    t(mat_top50)
  )
)

# Create treatment annotation

annotation_col = data.frame(
  Treatment = sample_info$condition
)

rownames(annotation_col) = rownames(sample_info)

# Make sure annotation follows the same sample order as the heatmap

annotation_col = annotation_col[
  colnames(mat_top50_scaled),
  ,
  drop = FALSE
]

# Treatment colors

annotation_colors = list(
  Treatment = c(
    Control = "#18C4CE",
    Cd = "#F7837D"
  )
)

# Save heatmap

png(
  file.path(
    plot_dir,
    "Top50_variable_genes_heatmap.png"
  ),
  width = 1200,
  height = 1500,
  res = 150
)

pheatmap(
  mat_top50_scaled,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 10,
  fontsize_col = 11,
  border_color = "grey80",
  main = "Top 50 Variable Genes",
  color = colorRampPalette(
    c(
      "#3B75AF",
      "#B9D8E8",
      "#FFFFD9",
      "#F9C477",
      "#F36F56"
    )
  )(100)
)

dev.off()


# Volcano plot

volcano_data = as.data.frame(res)

volcano_data$gene_id = rownames(
  volcano_data
)

volcano_data$significance = "Not significant"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange >= 1
] = "Upregulated"

volcano_data$significance[
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05 &
    volcano_data$log2FoldChange <= -1
] = "Downregulated"

volcano_data$neg_log10_padj = -log10(
  pmax(
    volcano_data$padj,
    .Machine$double.xmin
  )
)

volcano_plot = ggplot(
  volcano_data,
  aes(
    log2FoldChange,
    neg_log10_padj
  )
) +
  geom_point(
    aes(shape = significance),
    alpha = 0.7,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Volcano Plot - Cd vs Control",
    x = "log2 Fold Change",
    y = "-log10 adjusted P value"
  ) +
  theme_classic()

ggsave(
  file.path(
    plot_dir,
    "volcano_Cd_vs_Control.png"
  ),
  volcano_plot,
  width = 7,
  height = 6,
  dpi = 300
)


# Main enrichment analysis uses nuclear rice genes beginning with Os.
# ENSRNA features and chloroplast tRNA features are retained in the
# DEG tables but are not mixed into standard nuclear GO/KEGG enrichment.

universe_os = all_ids[grepl("^Os", all_ids)]
deg_os = rownames(sig)[grepl("^Os", rownames(sig))]
deg_up_os = rownames(deg_up)[grepl("^Os", rownames(deg_up))]
deg_down_os = rownames(deg_down)[grepl("^Os", rownames(deg_down))]

cat("\nNuclear Os genes\n")
cat("Background:", length(universe_os), "\n")
cat("DEGs:", length(deg_os), "\n")
cat("Up:", length(deg_up_os), "\n")
cat("Down:", length(deg_down_os), "\n")

# GO annotation from Ensembl Plants BioMart

# BioMart export:
# Oryza sativa Japonica Group genes (IRGSP-1.0)
# Selected attributes:
# Gene stable ID
# Gene name
# GO term accession
# GO term name
# GO term definition
# NCBI gene / Entrez ID

go_file = "C:/Users/nancy.AEC0/Downloads/mart_export (3).txt"

go_annot = read.delim(
  go_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

go_annot_clean = go_annot[
  !is.na(go_annot$GO.term.accession) &
    go_annot$GO.term.accession != "" &
    grepl("^Os", go_annot$Gene.stable.ID),
]

go_annot_clean = unique(
  go_annot_clean[
    ,
    c("Gene.stable.ID", "GO.term.accession", "GO.term.name")
  ]
)

term2gene = go_annot_clean[
  ,
  c("GO.term.accession", "Gene.stable.ID")
]

term2name = unique(
  go_annot_clean[
    ,
    c("GO.term.accession", "GO.term.name")
  ]
)

colnames(term2gene) = c("term", "gene")
colnames(term2name) = c("term", "name")

go_genes = unique(term2gene$gene)

universe_os_annot = intersect(universe_os, go_genes)
deg_os_annot = intersect(deg_os, go_genes)
deg_up_os_annot = intersect(deg_up_os, go_genes)
deg_down_os_annot = intersect(deg_down_os, go_genes)

cat("\nGO coverage\n")
cat("Background:", length(universe_os_annot), "\n")
cat("DEGs:", length(deg_os_annot), "\n")
cat("Up:", length(deg_up_os_annot), "\n")
cat("Down:", length(deg_down_os_annot), "\n")

# ------------------------------------------------------------
# 14. GO enrichment
# ------------------------------------------------------------

ego_all = enricher(
  gene = deg_os_annot,
  universe = universe_os_annot,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego_up = enricher(
  gene = deg_up_os_annot,
  universe = universe_os_annot,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego_down = enricher(
  gene = deg_down_os_annot,
  universe = universe_os_annot,
  TERM2GENE = term2gene,
  TERM2NAME = term2name,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ego_all_df = as.data.frame(ego_all)
ego_up_df = as.data.frame(ego_up)
ego_down_df = as.data.frame(ego_down)

write.csv(
  ego_all_df,
  file.path(result_dir, "GO_all_DEGs.csv"),
  row.names = FALSE
)
write.csv(
  ego_up_df,
  file.path(result_dir, "GO_upregulated_DEGs.csv"),
  row.names = FALSE
)
write.csv(
  ego_down_df,
  file.path(result_dir, "GO_downregulated_DEGs.csv"),
  row.names = FALSE
)

# GO plots
if (nrow(ego_all_df) > 0) {
  pdf(file.path(result_dir, "GO_all_DEGs_dotplot.pdf"), 9, 7)
  print(
    dotplot(
      ego_all,
      showCategory = min(20, nrow(ego_all_df)),
      font.size = 10
    ) +
      ggtitle("GO Enrichment - All DEGs")
  )
  dev.off()
}

if (nrow(ego_up_df) > 0) {
  pdf(file.path(result_dir, "GO_upregulated_dotplot.pdf"), 9, 7)
  print(
    dotplot(
      ego_up,
      showCategory = min(20, nrow(ego_up_df)),
      font.size = 10
    ) +
      ggtitle("GO Enrichment - Upregulated DEGs")
  )
  dev.off()
}

if (nrow(ego_down_df) > 0) {
  pdf(file.path(result_dir, "GO_downregulated_dotplot.pdf"), 9, 7)
  print(
    dotplot(
      ego_down,
      showCategory = min(20, nrow(ego_down_df)),
      font.size = 10
    ) +
      ggtitle("GO Enrichment - Downregulated DEGs")
  )
  dev.off()
}

# KEGG annotation

# KEGG organism code:
# dosa = Oryza sativa japonica (Japanese rice), RAPDB

kegg_link = read.delim(
  "https://rest.kegg.jp/link/pathway/dosa",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

kegg_list = read.delim(
  "https://rest.kegg.jp/list/pathway/dosa",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

term2gene_kegg = data.frame(
  term = sub("^dosa:", "", kegg_link$V2),
  gene = sub("^dosa:", "", kegg_link$V1),
  stringsAsFactors = FALSE
)

term2gene_kegg$term = gsub(
  "^path:",
  "",
  term2gene_kegg$term
)

term2name_kegg = data.frame(
  term = kegg_list$V1,
  name = sub(
    " - Oryza sativa japonica.*$",
    "",
    kegg_list$V2
  ),
  stringsAsFactors = FALSE
)

cat(
  "\nMatching KEGG pathways:",
  length(
    intersect(
      unique(term2gene_kegg$term),
      unique(term2name_kegg$term)
    )
  ),
  "\n"
)

# KEGG background

kegg_genes = unique(term2gene_kegg$gene)

universe_kegg = intersect(universe_os, kegg_genes)
deg_kegg = intersect(deg_os, kegg_genes)
deg_up_kegg = intersect(deg_up_os, kegg_genes)
deg_down_kegg = intersect(deg_down_os, kegg_genes)

term2gene_kegg_filtered = unique(
  term2gene_kegg[
    term2gene_kegg$gene %in% universe_kegg,
  ]
)

cat("\nKEGG coverage\n")
cat("Background:", length(universe_kegg), "\n")
cat("DEGs:", length(deg_kegg), "\n")
cat("Up:", length(deg_up_kegg), "\n")
cat("Down:", length(deg_down_kegg), "\n")

# KEGG enrichment
ekegg_all = enricher(
  gene = deg_kegg,
  universe = universe_kegg,
  TERM2GENE = term2gene_kegg_filtered,
  TERM2NAME = term2name_kegg,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ekegg_up = enricher(
  gene = deg_up_kegg,
  universe = universe_kegg,
  TERM2GENE = term2gene_kegg_filtered,
  TERM2NAME = term2name_kegg,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ekegg_down = enricher(
  gene = deg_down_kegg,
  universe = universe_kegg,
  TERM2GENE = term2gene_kegg_filtered,
  TERM2NAME = term2name_kegg,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ekegg_all_df = as.data.frame(ekegg_all)
ekegg_up_df = as.data.frame(ekegg_up)
ekegg_down_df = as.data.frame(ekegg_down)

write.csv(
  ekegg_all_df,
  file.path(result_dir, "KEGG_all_DEGs.csv"),
  row.names = FALSE
)
write.csv(
  ekegg_up_df,
  file.path(result_dir, "KEGG_upregulated_DEGs.csv"),
  row.names = FALSE
)
write.csv(
  ekegg_down_df,
  file.path(result_dir, "KEGG_downregulated_DEGs.csv"),
  row.names = FALSE
)

# KEGG plots
if (nrow(ekegg_up_df) > 0) {
  pdf(file.path(result_dir, "KEGG_upregulated_dotplot.pdf"), 9, 7)
  print(
    dotplot(
      ekegg_up,
      showCategory = min(20, nrow(ekegg_up_df)),
      font.size = 10
    ) +
      ggtitle("KEGG Enrichment - Upregulated DEGs")
  )
  dev.off()
}

if (nrow(ekegg_down_df) > 0) {
  pdf(file.path(result_dir, "KEGG_downregulated_dotplot.pdf"), 9, 7)
  print(
    dotplot(
      ekegg_down,
      showCategory = min(20, nrow(ekegg_down_df)),
      font.size = 10
    ) +
      ggtitle("KEGG Enrichment - Downregulated DEGs")
  )
  dev.off()
}

# Summary
summary_table = data.frame(
  Category = c(
    "DESeq2 features",
    "Significant DEGs",
    "Upregulated DEGs",
    "Downregulated DEGs",
    "Nuclear Os background",
    "Nuclear Os DEGs",
    "GO-annotated background",
    "GO-annotated DEGs",
    "KEGG background",
    "KEGG DEGs"
  ),
  Number = c(
    nrow(res),
    nrow(sig),
    nrow(deg_up),
    nrow(deg_down),
    length(universe_os),
    length(deg_os),
    length(universe_os_annot),
    length(deg_os_annot),
    length(universe_kegg),
    length(deg_kegg)
  )
)

write.csv(
  summary_table,
  file.path(result_dir, "analysis_summary.csv"),
  row.names = FALSE
)

print(summary_table)

cat("\nAnalysis completed successfully.\n")
cat("Results saved in:", result_dir, "\n")


print(summary_table)

cat("\nAnalysis completed successfully.\n")
cat("Results saved in:", result_dir, "\n")


cat(
  "\n============================================================\n"
)

cat(
  "RNA-seq DESeq2 + GO + KEGG pipeline finished successfully.\n"
)

cat(
  "============================================================\n"
)



