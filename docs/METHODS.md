# Recorded analysis methods

This description is based on the uploaded scripts and report. Final run logs take precedence where software versions or executed commands differ. The analysis was performed with WSL/Linux tools and a separate R workflow.

## Processing and alignment

Raw reads were assessed with FastQC and MultiQC. The archived Trimmomatic paired-end command uses `ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true`, `SLIDINGWINDOW:5:20` and `MINLEN:50`. The minimum length is 50 bases, not 150. The adapter file is the sequence collection selected for trimming; its name alone does not verify the original library preparation kit. The command does not explicitly implement a separate filter removing every read that contains an N. Trimmomatic's sliding-window behavior should not be described as an exact reproduction of a differently directed window algorithm in another method.

Retained paired reads were checked again and aligned with HISAT2 to an IRGSP-1.0 genome index. The archived HISAT2 command does not specify library strandedness or an annotation-derived splice-site file. SAMtools sorted and indexed the alignments and generated flagstat summaries. A missing flagstat report can be regenerated from an existing completed BAM without repeating alignment.

## Counting

The archived featureCounts command uses `-T 8 -p -B -C -t exon -g gene_id`. Although the Bash variable THREADS is 4, this specific command requests 8 threads. No explicit `-s` option is present, so the run uses the program's default unstranded setting. The appropriateness of that choice should be assessed using the library protocol or strandedness evidence.

The counting unit is version-dependent: from Subread 2.0.2 onward, `-p` declares paired-end input and `--countReadPairs` separately requests counting pairs. The archived command omits the latter. The exact version and featureCounts log must therefore be retained before calling these fragment counts. No counting-mode change was made during packaging.

## Differential expression

Count columns were matched by run accession to sample IDs and reordered to metadata. Features were retained with at least 10 counts in at least three libraries. DESeq2 used `~ condition`, with Control as the reference and the contrast Cd versus Control. Significance was defined as BH-adjusted P < 0.05 and absolute log₂ fold change ≥ 1.

The archived R script calls `results()` without an explicit `alpha`, while selecting DEGs afterward at padj < 0.05. Record this as executed rather than silently changing the call; alpha can influence independent filtering. It also rounds the input counts. featureCounts inputs should already be nonnegative integer counts; rounding must not be used to substitute TPM or other normalized expression values for counts.

Variance-stabilized values (`vst`, blind = FALSE) were used for QC. `plotPCA` selects its default most-variable features (normally up to 500), sample correlations use Pearson correlation, and sample distances use Euclidean distance. The top-50 heatmap selects features by variance and scales each row. These features are not necessarily statistically significant DEGs.

The archived distance heatmap calls pheatmap on a distance matrix without explicitly passing the original distance object for clustering. Its tree therefore reflects the plotting function's clustering of that supplied matrix. Interpret the displayed pairwise distances directly; document any later change to the clustering arguments.

## Enrichment

GO and KEGG analyses were conducted for all, upregulated and downregulated significant nuclear `Os` identifiers using clusterProfiler's `enricher`. GO mappings came from Ensembl Plants BioMart; KEGG mappings were requested using organism code `dosa`. BH correction, pvalueCutoff = 0.05 and qvalueCutoff = 0.05 are specified.

The archived background consists of `Os` identifiers in the DESeq2 result object intersected with each annotation resource. It does not explicitly restrict the background to features with non-missing adjusted P values. Describe this accurately; a revision to the eligible background requires rerunning enrichment and replacing its outputs. Archive annotation retrieval dates and exact input mappings where redistribution is permitted. Live database queries can change over time.
