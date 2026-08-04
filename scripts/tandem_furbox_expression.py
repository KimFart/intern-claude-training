#!/usr/bin/env python3
"""Build the gene-level log2FC table for the tandem Fur-box hypothesis (H2).

For each Fur-regulated gene (from site_to_gene.py's mapping), sums per-base
RNA-seq depth over the gene body from the four *_linear.gff tracks (built by
notebooks/makegff.py WITHOUT --log_scale, so depth is linear and summable --
see the plan's note on why the existing --log_scale GFFs can't be used for
this), normalizes each replicate by its own total mapped-read count, averages
the two replicates per condition, and computes log2FC(starve / replete).

Each *_linear.gff row is a single non-zero-depth base position (start == end,
one row per base -- notebooks/makegff.py does not run-length-encode), so a
gene-body sum is just a strand-matched boolean mask + sum over the score column.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

DATA_REF = "/workspaces/intern-claude-training/data/reference"
DATA_MOD5 = "/workspaces/intern-claude-training/data/module5"
SITE_TO_GENE_PATH = os.path.join(DATA_MOD5, "site_to_gene.tsv")
OUT_PATH = os.path.join(DATA_MOD5, "tandem_furbox_hypothesis_table.tsv")

# SRR -> (condition, replicate)
REPLICATES = {
    "SRR1168133": ("replete", 1),
    "SRR1168134": ("replete", 2),
    "SRR1168135": ("starve", 1),
    "SRR1168136": ("starve", 2),
}
PSEUDOCOUNT = 1.0
CPM_SCALE = 1e6


def total_mapped_reads(srr):
    """Total mapped reads for a sorted BAM, via samtools idxstats (sum of column 3)."""
    bam_path = os.path.join(DATA_REF, f"{srr}_sorted.bam")
    if not os.path.isfile(bam_path):
        sys.exit(f"Error: sorted BAM not found: {bam_path}")
    result = subprocess.run(
        ["samtools", "idxstats", bam_path], capture_output=True, text=True, check=True
    )
    total = sum(int(line.split("\t")[2]) for line in result.stdout.strip().splitlines())
    return total


def load_linear_gff(srr):
    """Load a *_linear.gff (per-base depth, one row per non-zero position)."""
    gff_path = os.path.join(DATA_REF, f"{srr}_linear.gff")
    if not os.path.isfile(gff_path):
        sys.exit(
            f"Error: {gff_path} not found. Run notebooks/makegff.py --flip "
            f"--separate_strand (no --log_scale) on {srr}_sorted.bam first."
        )
    df = pd.read_csv(
        gff_path, sep="\t", header=None, usecols=[3, 4, 5, 6],
        names=["start", "end", "score", "strand"],
        dtype={"start": "int32", "end": "int32", "score": "float32", "strand": "category"},
    )
    return df


def gene_depth_sums(gff_df, genes):
    """Sum of depth per gene (matching strand, position within [start, end])."""
    sums = {}
    for strand in ("+", "-"):
        strand_rows = gff_df[gff_df["strand"] == strand]
        positions = strand_rows["start"].to_numpy()
        scores = strand_rows["score"].to_numpy()
        for _, gene in genes[genes["matched_strand"] == strand].iterrows():
            mask = (positions >= gene["matched_start"]) & (positions <= gene["matched_end"])
            sums[gene["matched_locus_tag"]] = float(scores[mask].sum())
    return sums


def main():
    if not os.path.isfile(SITE_TO_GENE_PATH):
        sys.exit(f"Error: {SITE_TO_GENE_PATH} not found. Run scripts/site_to_gene.py first.")
    sites = pd.read_csv(SITE_TO_GENE_PATH, sep="\t")

    genes = sites.drop_duplicates(subset=["matched_locus_tag"])[
        ["matched_locus_tag", "matched_gene", "matched_start", "matched_end", "matched_strand"]
    ]
    print(f"Unique genes to quantify: {len(genes)} (from {len(sites)} peaks)")

    normalized = {}  # locus_tag -> {"replete": [rep1, rep2], "starve": [rep1, rep2]}
    for srr, (condition, rep) in REPLICATES.items():
        print(f"Processing {srr} ({condition} rep{rep})...")
        mapped_reads = total_mapped_reads(srr)
        gff_df = load_linear_gff(srr)
        depth_sums = gene_depth_sums(gff_df, genes)
        for locus_tag, depth_sum in depth_sums.items():
            cpm = depth_sum * CPM_SCALE / mapped_reads
            normalized.setdefault(locus_tag, {"replete": [], "starve": []})[condition].append(cpm)
        del gff_df

    rows = []
    for _, gene in genes.iterrows():
        locus_tag = gene["matched_locus_tag"]
        replete_reps = normalized[locus_tag]["replete"]
        starve_reps = normalized[locus_tag]["starve"]
        replete_mean = float(np.mean(replete_reps))
        starve_mean = float(np.mean(starve_reps))
        log2fc = np.log2((starve_mean + PSEUDOCOUNT) / (replete_mean + PSEUDOCOUNT))
        rows.append({
            "matched_locus_tag": locus_tag,
            "matched_gene": gene["matched_gene"],
            "replete_rep1_cpm": replete_reps[0],
            "replete_rep2_cpm": replete_reps[1],
            "starve_rep1_cpm": starve_reps[0],
            "starve_rep2_cpm": starve_reps[1],
            "replete_mean_cpm": replete_mean,
            "starve_mean_cpm": starve_mean,
            "log2FC": log2fc,
        })
    gene_table = pd.DataFrame(rows)

    combined = sites.merge(gene_table, on=["matched_locus_tag", "matched_gene"], how="left")
    os.makedirs(DATA_MOD5, exist_ok=True)
    combined.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"Hypothesis table written to: {OUT_PATH} ({len(combined)} rows, {len(gene_table)} unique genes)")


if __name__ == "__main__":
    main()
