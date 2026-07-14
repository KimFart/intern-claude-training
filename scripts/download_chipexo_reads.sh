#!/usr/bin/env bash
# Exercise 4 — Download a single-end ChIP-exo sample from GEO series GSE54901
# (Seo et al. 2014, Fur ChIP-exo in E. coli K-12 MG1655) into data/reference/.
#
# GSE54901 is a SuperSeries with 12 ChIP-exo samples + 8 RNA-seq samples.
# You need one of the ChIP-exo (not RNA-seq) SRR accessions — look one up via
# the SRA Run Selector: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE54899
# (the ChIP-exo SubSeries — the SuperSeries GSE54901 itself isn't indexed there)
# then pass it as the first argument to this script, e.g.:
#   ./scripts/download_chipexo_reads.sh SRRxxxxxxx [threads]
#
# Uses prefetch + fasterq-dump (not plain fastq-dump): a real ChIP-exo run is
# 10-50 million reads, and fastq-dump fetches/converts single-threaded, row by
# row, over the network. prefetch downloads the .sra file once (resumable if
# the connection drops); fasterq-dump then converts it locally with multiple
# threads, which is far faster at that scale.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 <SRR_accession> [threads]" >&2
    echo "  SRR_accession   an SRA run accession for one of GSE54901's ChIP-exo samples" >&2
    echo "  threads         thread count for fasterq-dump (default: nproc)" >&2
    exit 1
fi

SRR="$1"
THREADS="${2:-$(nproc)}"
OUTPUT_DIR="data/reference"
SRA_CACHE_DIR="$OUTPUT_DIR/.sra_cache"
PREFETCH="/opt/sratoolkit/bin/prefetch"
FASTERQ_DUMP="/opt/sratoolkit/bin/fasterq-dump"

mkdir -p "$OUTPUT_DIR" "$SRA_CACHE_DIR"

echo "Prefetching $SRR into $SRA_CACHE_DIR"

# --output-directory   where prefetch stores <SRR>/<SRR>.sra; kept separate
#                       from data/reference/ so the .sra cache doesn't get
#                       mistaken for pipeline input
"$PREFETCH" "$SRR" --output-directory "$SRA_CACHE_DIR"

echo "Converting $SRR to FASTQ with $THREADS threads"

# --outdir    write the .fastq file into data/reference/
# --threads   parallelize the SRA -> FASTQ conversion (fastq-dump has no
#             equivalent flag; it's single-threaded, which is the whole
#             reason for switching to fasterq-dump for larger runs)
# --progress  show a progress bar - useful for multi-million-read runs
# (no --split-files: this is single-end data, so there's only one read per
#  spot; --split-files only matters for paired-end runs with _1/_2 mates)
"$FASTERQ_DUMP" "$SRA_CACHE_DIR/$SRR/$SRR.sra" \
    --outdir "$OUTPUT_DIR" \
    --threads "$THREADS" \
    --progress

echo "Cleaning up .sra cache for $SRR"
rm -rf "${SRA_CACHE_DIR:?}/$SRR"

echo "Done: $OUTPUT_DIR/${SRR}.fastq"
