#!/usr/bin/env bash
set -Eeuo pipefail
# =============================================================================
# RNA-seq ANALYSIS PIPELINE
#
# Study:
#   Transcriptome analysis of rice (Oryza sativa L.) shoots responsive to
#   cadmium stress
#
# Dataset:
#   PRJNA544413
#
# Experimental design:
#   Control: 0 µM Cd, 7 days      (3 biological replicates)
#   Cd stress: 75 µM Cd, 7 days   (3 biological replicates)
#   Paired-end RNA-seq
#
# Workflow:
#   1.  FastQC                  - quality control of raw reads
#   2.  MultiQC                 - combined raw-read report
#   3.  Trimmomatic             - adapter and quality trimming
#   4.  FastQC                  - quality control of trimmed reads
#   5.  MultiQC                 - combined trimmed-read report
#   6.  HISAT2                  - prepare rice reference genome
#   7.  HISAT2                  - align reads to the rice genome
#   8.  SAMtools                - convert and sort SAM files
#   9.  SAMtools                - index BAM files
#  10.  SAMtools                - alignment quality statistics
#  11.  featureCounts           - gene-level read counting
#
# The final output of this Bash script is:
#   05_featurecounts/gene_counts.txt
#
# All commands are run under WSL/Linux.
# =============================================================================


# =============================================================================
# PROJECT SETUP
# =============================================================================

PROJECT="/mnt/d/Nancy/Project"
THREADS=4

# Activate the Conda environment containing the RNA-seq tools.
conda activate rnaseq

# Move into the project directory.
cd "$PROJECT" || exit 1


# CREATE PROJECT DIRECTORIES

echo ""
echo "Preparing project directories"
mkdir -p 01_fastqc_raw
mkdir -p 02_multiqc_raw
mkdir -p 03_trimmomatic
mkdir -p 04_fastqc_trimmed
mkdir -p 05_multiqc_trimmed
mkdir -p 04_hisat2/sam
mkdir -p 04_hisat2/bam
mkdir -p 05_featurecounts
mkdir -p 06_logs

echo "Project directories are ready."

# CHECK Raw FASTQ Files
# =============================================================================

echo ""
echo "Checking raw FASTQ files"
echo "Raw FASTQ files:"
ls -lh *.fastq

echo ""
echo "Number of FASTQ files:"
ls *_1.fastq *_2.fastq 2>/dev/null | wc -l

echo ""
echo "Expected samples:"
echo "  SRR9113519  Control"
echo "  SRR9113520  Control"
echo "  SRR9113521  Control"
echo "  SRR9113522  Cd"
echo "  SRR9113523  Cd"
echo "  SRR9113524  Cd"
echo ""

# Check that every expected paired-end file is present.
SAMPLES=(
    SRR9113519
    SRR9113520
    SRR9113521
    SRR9113522
    SRR9113523
    SRR9113524
)

for SAMPLE in "${SAMPLES[@]}"; do

    if [ ! -f "${SAMPLE}_1.fastq" ]; then
        echo "ERROR: Missing ${SAMPLE}_1.fastq"
        exit 1
    fi

    if [ ! -f "${SAMPLE}_2.fastq" ]; then
        echo "ERROR: Missing ${SAMPLE}_2.fastq"
        exit 1
    fi

done
echo "All expected FASTQ files were found."


# QC Of Raw Reads
# =============================================================================

echo ""
echo "FastQC - raw reads"

fastqc \
    -t "$THREADS" \
    -o 01_fastqc_raw \
    *.fastq

echo ""
echo "Raw-read FastQC completed."


# COMBINE RAW FASTQC RESULTS WITH MULTIQC
# =============================================================================

echo ""
echo "MultiQC - raw reads"

multiqc \
    01_fastqc_raw \
    -o 02_multiqc_raw

echo ""
echo "Raw-read MultiQC report completed."


# ADAPTER AND QUALITY TRIMMING
# =============================================================================

echo ""
echo "Trimmomatic preprocessing"

# Find the TruSeq paired-end adapter file inside the active Conda environment.
ADAPTER=$(find "$CONDA_PREFIX" \
    -name "TruSeq3-PE.fa" \
    | head -n 1)

# Stop if the adapter file cannot be found.
if [ -z "$ADAPTER" ]; then
    echo "ERROR: TruSeq3-PE.fa was not found."
    echo "Please check your Trimmomatic installation."
    exit 1
fi

echo "Adapter file:"
echo "  $ADAPTER"
echo ""

# Process each paired-end sample.
# =============================================================================

for R1 in *_1.fastq; do

    SAMPLE="${R1%_1.fastq}"
    R2="${SAMPLE}_2.fastq"
    echo "Trimming sample: $SAMPLE"

    # Make sure the reverse-read file exists.
    if [ ! -f "$R2" ]; then
        echo "ERROR: Missing paired-end file: $R2"
        exit 1
    fi

    # Trimmomatic settings:
    # ILLUMINACLIP:
    #   Remove Illumina adapter sequences.
    # SLIDINGWINDOW:5:20
    #   Trim when the average quality in a 5-base window falls below Q20.
    # MINLEN:50
    #   Discard reads shorter than 50 bp after trimming.

    trimmomatic PE \
        -threads "$THREADS" \
        -phred33 \
        "$R1" \
        "$R2" \
        "03_trimmomatic/${SAMPLE}_1_paired.fastq" \
        "03_trimmomatic/${SAMPLE}_1_unpaired.fastq" \
        "03_trimmomatic/${SAMPLE}_2_paired.fastq" \
        "03_trimmomatic/${SAMPLE}_2_unpaired.fastq" \
        ILLUMINACLIP:"$ADAPTER":2:30:10:8:true \
        SLIDINGWINDOW:5:20 \
        MINLEN:50 \
        > "06_logs/${SAMPLE}_trimmomatic.log" 2>&1

    echo "Finished trimming: $SAMPLE"
    echo ""

done

echo "Trimmomatic preprocessing completed."


# QC of Trimmed Reads
# =============================================================================

echo ""
echo "FastQC - trimmed reads"

fastqc \
    -t "$THREADS" \
    -o 04_fastqc_trimmed \
    03_trimmomatic/*_paired.fastq

echo ""
echo "Trimmed-read FastQC completed."


# MULTIQC Report for Trimmed Reads
# =============================================================================

echo ""
echo "============================================================"
echo "MultiQC - trimmed reads"
echo "============================================================"

multiqc \
    04_fastqc_trimmed \
    -o 05_multiqc_trimmed

echo ""
echo "Trimmed-read MultiQC report completed."


# PREPARE THE RICE REFERENCE GENOME
# =============================================================================

echo ""
echo "============================================================"
echo "Preparing the rice reference genome"
echo "============================================================"

cd "$PROJECT/annotation" || exit 1

# The rice reference genome is provided as a BGZF-compressed FASTA file.
# Convert it to a standard FASTA file for HISAT2 indexing.

if [ ! -f "unmasked.fa.bgz" ]; then
    echo "ERROR: Reference genome unmasked.fa.bgz was not found."
    exit 1
fi

if [ ! -f "IRGSP-1.0.fa" ]; then

    echo "Creating IRGSP-1.0.fa from unmasked.fa.bgz..."

    gzip -dc unmasked.fa.bgz > IRGSP-1.0.fa

else

    echo "IRGSP-1.0.fa already exists. Using the existing FASTA file."

fi

echo ""
echo "Reference genome:"
ls -lh IRGSP-1.0.fa


# BUILD THE HISAT2 REFERENCE INDEX
# =============================================================================

echo ""
echo "Building HISAT2 genome index"
echo "============================================================"

# Build the index only if it does not already exist.
if [ ! -f "IRGSP-1.0.1.ht2" ]; then

    hisat2-build \
        -p "$THREADS" \
        IRGSP-1.0.fa \
        IRGSP-1.0

    echo ""
    echo "HISAT2 index successfully created."

else

    echo ""
    echo "HISAT2 index already exists. Skipping index construction."

fi

echo ""
echo "HISAT2 index files:"
ls -lh IRGSP-1.0*.ht2


# Align Trimmed Reads to the Rice Genome Using HISAT2
# =============================================================================

echo ""
echo "HISAT2 alignment"
echo "============================================================"

cd "$PROJECT" || exit 1

for SAMPLE in "${SAMPLES[@]}"; do

    echo ""
    echo "Aligning sample: $SAMPLE"
    echo "------------------------------------------------------------"

    hisat2 \
        -p "$THREADS" \
        -x "$PROJECT/annotation/IRGSP-1.0" \
        -1 "$PROJECT/03_trimmomatic/${SAMPLE}_1_paired.fastq" \
        -2 "$PROJECT/03_trimmomatic/${SAMPLE}_2_paired.fastq" \
        -S "$PROJECT/04_hisat2/sam/${SAMPLE}.sam" \
        > "$PROJECT/06_logs/${SAMPLE}_hisat2.log" 2>&1

    echo "Alignment completed: $SAMPLE"

done

# Convert & Sort SAM Files
# =============================================================================

echo ""
echo "============================================================"
echo "Converting and sorting SAM files"
echo "============================================================"

cd "$PROJECT/04_hisat2/sam" || exit 1

for SAMPLE in "${SAMPLES[@]}"; do

    echo "Sorting: $SAMPLE"

    samtools sort \
        -@ "$THREADS" \
        -o "$PROJECT/04_hisat2/bam/${SAMPLE}.sorted.bam" \
        "${SAMPLE}.sam"

    echo "Finished sorting: $SAMPLE"

done


# Index sorted BAM files
# =============================================================================

echo ""
echo "Indexing BAM files"
echo "============================================================"

cd "$PROJECT/04_hisat2/bam" || exit 1

for SAMPLE in "${SAMPLES[@]}"; do

    echo "Indexing: $SAMPLE"

    samtools index \
        -@ "$THREADS" \
        "${SAMPLE}.sorted.bam"

    echo "Finished indexing: $SAMPLE"

done

# Alignment Quality Assessment
# =============================================================================

echo ""
echo "Alignment statistics"
echo "============================================================"

cd "$PROJECT/04_hisat2/bam" || exit 1

for SAMPLE in "${SAMPLES[@]}"; do

    echo ""
    echo "Alignment statistics: $SAMPLE"
    echo "------------------------------------------------------------"

    samtools flagstat \
        "${SAMPLE}.sorted.bam" \
        > "$PROJECT/06_logs/${SAMPLE}_flagstat.txt"

    cat "$PROJECT/06_logs/${SAMPLE}_flagstat.txt"

done

# GENE-LEVEL READ COUNTING
# =============================================================================

echo ""
echo "FeatureCounts - gene-level read counting"
echo "============================================================"

# Check that the annotation file exists before counting.
if [ ! -f "$PROJECT/annotation/genes.gtf" ]; then
    echo "ERROR: GTF annotation file was not found:"
    echo "  $PROJECT/annotation/genes.gtf"
    exit 1
fi

# Count reads assigned to exons and summarize them at the gene level.
#
# -T "$THREADS"
#   Use $THREADS CPU threads for featureCounts.
#
# -p
#   Paired-end reads.
#
# -B
#   Require both ends of a pair to be successfully aligned.
#
# -C
#   Exclude chimeric fragments.
#
# -t exon
#   Count exon features.
#
# -g gene_id
#   Summarize exon counts by gene_id.
#
# -a
#   Rice GTF annotation.
#
# -o
#   Final gene-level count matrix.

featureCounts \
    -T "$THREADS" \
    -p \
    -B \
    -C \
    -t exon \
    -g gene_id \
    -a "$PROJECT/annotation/genes.gtf" \
    -o "$PROJECT/05_featurecounts/gene_counts.txt" \
    "$PROJECT/04_hisat2/bam/SRR9113519.sorted.bam" \
    "$PROJECT/04_hisat2/bam/SRR9113520.sorted.bam" \
    "$PROJECT/04_hisat2/bam/SRR9113521.sorted.bam" \
    "$PROJECT/04_hisat2/bam/SRR9113522.sorted.bam" \
    "$PROJECT/04_hisat2/bam/SRR9113523.sorted.bam" \
    "$PROJECT/04_hisat2/bam/SRR9113524.sorted.bam" \
    > "$PROJECT/06_logs/featureCounts.log" 2>&1

echo ""
echo "featureCounts completed."


# Check the Final Count Matrix
# =============================================================================

echo ""
echo "Checking featureCounts output"
echo "============================================================"

if [ ! -f "$PROJECT/05_featurecounts/gene_counts.txt" ]; then
    echo "ERROR: featureCounts did not produce the expected count file."
    exit 1
fi

echo ""
echo "Final gene-level count matrix:"
ls -lh "$PROJECT/05_featurecounts/gene_counts.txt"

echo ""
echo "First five lines:"
head -n 5 "$PROJECT/05_featurecounts/gene_counts.txt"


# Pipeline COMPLETED
# =============================================================================

echo ""
echo "============================================================"
echo "RNA-seq preprocessing and gene counting completed"
echo "============================================================"
echo ""
echo "Final featureCounts output:"
echo "  $PROJECT/05_featurecounts/gene_counts.txt"
