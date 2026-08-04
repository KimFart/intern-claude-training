#!/usr/bin/env bash
# One-shot reproduction of the Module 5 tandem Fur-box hypothesis test.
#
# 1. Download + align + sort + index SRR1168134 and SRR1168136 (the second
#    biological replicate of each RNA-seq condition), reusing the same
#    paired-end recipe Module 4 used for SRR1168133/SRR1168135.
# 2. Regenerate all four RNA-seq GFFs with notebooks/makegff.py (the
#    lab-supplied script, NOT scripts/makegff.py) WITHOUT --log_scale, since
#    the hypothesis-table step needs linear, summable depth.
# 3. Run the Python analysis: site-to-gene mapping, gene-level log2FC table,
#    and the H1/H2 statistical tests + figures.
#
# Every step is skipped if its output already exists, so re-running after a
# partial failure just picks up where it left off.

set -euo pipefail
cd "$(dirname "$0")/.."

OUT="data/reference"
THREADS="${THREADS:-$(nproc)}"
export PATH="/opt/sratoolkit/bin:$PATH"

fetch_replicate() {
    local srr="$1"
    if [[ -f "$OUT/${srr}_sorted.bam" && -f "$OUT/${srr}_sorted.bam.bai" ]]; then
        echo "[$srr] Sorted BAM + index already exist, skipping."
        return
    fi
    if [[ ! -f "$OUT/${srr}_1.fastq" || ! -f "$OUT/${srr}_2.fastq" ]]; then
        echo "[$srr] Downloading paired-end reads..."
        fasterq-dump --split-files "$srr" -O "$OUT/" -e "$THREADS"
    fi
    echo "[$srr] Aligning (paired-end, NC_000913.2_index)..."
    bowtie2 --very-fast -X 1000 -3 3 -p "$THREADS" --no-mixed --no-discordant \
        -x "$OUT/NC_000913.2_index" \
        -1 "$OUT/${srr}_1.fastq" -2 "$OUT/${srr}_2.fastq" \
        -S "$OUT/${srr}.sam"
    samtools view -bS "$OUT/${srr}.sam" -o "$OUT/${srr}.bam"
    samtools sort -@ "$THREADS" "$OUT/${srr}.bam" -o "$OUT/${srr}_sorted.bam"
    samtools index "$OUT/${srr}_sorted.bam"
    echo "[$srr] Done."
}

make_linear_gff() {
    local srr="$1"
    local out_gff="$OUT/${srr}_linear.gff"
    if [[ -f "$out_gff" ]]; then
        echo "[$srr] $out_gff already exists, skipping."
        return
    fi
    if [[ ! -f "$OUT/${srr}_sorted.bam" ]]; then
        echo "Error: $OUT/${srr}_sorted.bam not found -- run fetch_replicate first." >&2
        exit 1
    fi
    echo "[$srr] Generating linear (non-log-scaled) RNA-seq GFF..."
    python notebooks/makegff.py --flip --separate_strand "$OUT/${srr}_sorted.bam" "$out_gff"
}

echo "=== Step 1: acquire the second replicate ==="
fetch_replicate SRR1168134
fetch_replicate SRR1168136

echo "=== Step 2: regenerate linear GFFs for all 4 RNA-seq samples ==="
for srr in SRR1168133 SRR1168134 SRR1168135 SRR1168136; do
    if [[ ! -f "$OUT/${srr}_sorted.bam" ]]; then
        echo "Error: $OUT/${srr}_sorted.bam not found. SRR1168133/SRR1168135 are expected " \
             "to already exist from Module 4." >&2
        exit 1
    fi
    make_linear_gff "$srr"
done

echo "=== Step 3: site-to-gene mapping ==="
python scripts/site_to_gene.py

echo "=== Step 4: gene-level log2FC table ==="
python scripts/tandem_furbox_expression.py

echo "=== Step 5: H1/H2 statistical tests + figures ==="
(cd scripts && python tandem_furbox_hypothesis_test.py)

echo "Pipeline complete. Outputs in data/module5/:"
ls -la data/module5/
