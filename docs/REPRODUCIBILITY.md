# Reproducibility status

## Available in this draft

- Methods and results descriptions grounded in the uploaded report and scripts.
- Ten result images extracted from the report, with provenance.
- Exact copies of the older uploaded Bash and R scripts under `scripts/archive/`.
- A documented run-to-sample mapping from the user's supplied metadata text.

## Still needed from the final analysis

1. The actual metadata TSV and the exact scripts used for the final run.
2. The unnormalized gene count matrix, featureCounts assignment summary and counting log.
3. Full DESeq2 and enrichment CSV outputs, not only screenshots or selected significant genes.
4. Software versions, R session information, reference/annotation sources and retrieval dates.
5. The final QC reports and alignment logs for all six libraries.
6. Scripts used to create any publication-specific figures absent from the archived R file, including the fold-enrichment bar charts and annotated heatmaps.

The repository is not yet a verified reproduction environment. No raw-read alignment, counting, DESeq2 or enrichment rerun was performed during packaging. The final outputs cannot be reconstructed from figure pixels.

## Resolve before describing an exact reproduction

- The old Bash `echo` labels disagree with the later metadata. Those messages do not set the DESeq2 condition, but the actual metadata must be confirmed.
- The archived R file points to `05_featurecounts/sample_info.tsv`, while the user reported using `D:/Nancy/Project/metadata.tsv`.
- The archived BioMart input path points to a personal Downloads folder. Record the actual file and replace the path in a documented publication copy.
- featureCounts counting mode depends on the version and the presence of `--countReadPairs`; strandedness was not explicitly set.
- Exact reference files and software versions were not supplied.
- The archived R script queries KEGG live. Preserve the retrieval date and input mapping provenance.
- Some report figures include formatting or calculations not present in the archived R file. Retain the additional plotting code or identify those figures as report-derived illustrations.

## Minimal route to a reusable release

First preserve the exact executed scripts and outputs. Then make documented path-only edits to a separate publication copy. Readers may start from the published count matrix for differential expression once the metadata and R environment are supplied. Rebuilding the read-processing steps also requires the public reads, matching reference and annotation, and the recorded command versions/settings.

If a methodological correction is made, record it as a new analysis version and replace the affected downstream tables and figures. Do not present old results as outputs of newly changed commands. The original computational settings can be documented honestly without claiming that every possible alternative would give the same result.
