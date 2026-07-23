#!/usr/bin/env bash
# Exercise 4b — Download a reference genome FASTA by NCBI accession into
# data/reference/<accession>.fasta, using entrez-direct's efetch.
#
# This is what bowtie2-build (Exercise 7, step 1) turns into the alignment
# index, so the accession must match the genome the reads were aligned from
# (e.g. NC_000913.3 for E. coli K-12 MG1655).

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <accession>" >&2
    echo "  accession   an NCBI nucleotide accession, e.g. NC_000913.3" >&2
    exit 1
fi

ACCESSION="$1"
OUTPUT_DIR="data/reference"
OUTPUT_FASTA="$OUTPUT_DIR/${ACCESSION}.fasta"
EFETCH="/opt/conda/envs/sbml/bin/efetch"

mkdir -p "$OUTPUT_DIR"

echo "Downloading $ACCESSION as FASTA into $OUTPUT_FASTA"

# -db nucleotide   query NCBI's nucleotide database
# -id              the accession to fetch (e.g. NC_000913.3)
# -format fasta    return the record as FASTA (header + sequence), not a GenBank flat file
"$EFETCH" -db nucleotide -id "$ACCESSION" -format fasta > "$OUTPUT_FASTA"

echo "Done: $OUTPUT_FASTA"
