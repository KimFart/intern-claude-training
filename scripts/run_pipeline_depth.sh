#!/usr/bin/env bash
# ChIP-exo single-end pipeline: FASTQ + reference FASTA -> sorted, indexed BAM -> GFF.
# Same as run_pipeline.sh, but the final GFF is scored with makegff_depth.py's
# 5'-position read depth instead of a hardcoded '.'.
# <reference>/<reads> may be a local file path or an accession to download if missing.

set -euo pipefail

usage() {
    echo "Usage: $0 <reference> <reads> <sample_name> [output_dir] [threads]" >&2
    echo "  reference   path to a reference FASTA, or an NCBI accession (e.g. NC_000913.3) to fetch" >&2
    echo "  reads       path to single-end FASTQ reads, or an SRA accession (e.g. SRR1177157) to fetch" >&2
    echo "  sample_name sample identifier, used to name all output files" >&2
    echo "  output_dir  directory for downloaded/generated files (default: data/reference)" >&2
    echo "  threads     thread count for bowtie2/samtools sort (default: nproc)" >&2
    exit 1
}

if [[ $# -lt 3 || $# -gt 5 ]]; then
    usage
fi

REFERENCE_ARG="$1"
READS_ARG="$2"
SAMPLE="$3"
OUTPUT_DIR="${4:-data/reference}"
THREADS="${5:-$(nproc)}"

mkdir -p "$OUTPUT_DIR"

resolve_fastq_dump() {
    if command -v fastq-dump &>/dev/null; then
        echo "fastq-dump"
    elif [[ -x /opt/sratoolkit/bin/fastq-dump ]]; then
        echo "/opt/sratoolkit/bin/fastq-dump"
    else
        echo "Error: fastq-dump not found on PATH or at /opt/sratoolkit/bin/fastq-dump" >&2
        exit 1
    fi
}

# --- Resolve reference: existing path, previously-downloaded file, or fetch by accession ---
if [[ -f "$REFERENCE_ARG" ]]; then
    REFERENCE_FASTA="$REFERENCE_ARG"
else
    REFERENCE_FASTA="$OUTPUT_DIR/${REFERENCE_ARG}.fasta"
    if [[ -f "$REFERENCE_FASTA" ]]; then
        echo "Reference already downloaded: $REFERENCE_FASTA"
    else
        echo "Downloading reference genome $REFERENCE_ARG from NCBI"
        efetch -db nucleotide -id "$REFERENCE_ARG" -format fasta > "$REFERENCE_FASTA"
    fi
fi

# --- Resolve reads: existing path, previously-downloaded file, or fetch by accession ---
if [[ -f "$READS_ARG" ]]; then
    FASTQ="$READS_ARG"
else
    FASTQ="$OUTPUT_DIR/${READS_ARG}.fastq"
    if [[ -f "$FASTQ" ]]; then
        echo "Reads already downloaded: $FASTQ"
    else
        echo "Downloading reads $READS_ARG from SRA"
        FASTQ_DUMP="$(resolve_fastq_dump)"
        "$FASTQ_DUMP" "$READS_ARG" --outdir "$OUTPUT_DIR"
    fi
fi

INDEX_BASE="$OUTPUT_DIR/${SAMPLE}_index"
SAM="$OUTPUT_DIR/${SAMPLE}.sam"
BAM="$OUTPUT_DIR/${SAMPLE}.bam"
SORTED_BAM="$OUTPUT_DIR/${SAMPLE}_sorted.bam"
GFF="$OUTPUT_DIR/${SAMPLE}_chipexo_depth.gff"

# --- Step 1: bowtie2-build (skip if index already exists) ---
if [[ -f "${INDEX_BASE}.1.bt2" || -f "${INDEX_BASE}.1.bt2l" ]]; then
    echo "Index already exists at ${INDEX_BASE}.*, skipping bowtie2-build."
else
    echo "Building bowtie2 index: $INDEX_BASE"
    bowtie2-build "$REFERENCE_FASTA" "$INDEX_BASE"
fi

# --- Step 2: bowtie2 align (single-end) ---
echo "Aligning $FASTQ to $INDEX_BASE using $THREADS threads"
bowtie2 -p "$THREADS" -x "$INDEX_BASE" -U "$FASTQ" -S "$SAM"

# --- Step 3: SAM -> BAM ---
echo "Converting SAM to BAM"
samtools view -bS "$SAM" -o "$BAM"

# --- Step 4: sort ---
# samtools sort's -@ is *additional* threads on top of the main one (unlike
# bowtie2's -p, which is the total), so this uses one more thread than bowtie2 did.
echo "Sorting BAM using $THREADS additional threads"
samtools sort -@ "$THREADS" "$BAM" -o "$SORTED_BAM"

# --- Step 5: index ---
echo "Indexing sorted BAM"
samtools index "$SORTED_BAM"

# --- Step 6: BAM -> GFF, scored by 5'-position read depth ---
echo "Generating depth-scored GFF: $GFF"
python scripts/makegff_depth.py "$SORTED_BAM" "$GFF"

echo "Pipeline complete."
echo "  Sorted BAM: $SORTED_BAM"
echo "  Index:      ${SORTED_BAM}.bai"
echo "  GFF:        $GFF"
