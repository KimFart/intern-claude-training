#!/usr/bin/env python3
"""Map each regulatory-location Fur binding site to its nearest annotated gene.

Data sources:
- Binding sites: ncomms5910_supplementary_data1.nomerge.xlsx (Seo et al. 2014
  Supplementary Data 1, merged-cell-fixed so every row has its own
  Transcription Unit / Location / Mode / '# of motifs' value)
- Gene coordinates: ec_annotation_20100903_DHK_cSRNA_with_ortho.gff

Transcription Unit names in the table are operon-letter-suffix strings
(e.g. 'fhuACDB', 'entCEBAH') that are ambiguous to split into individual
gene names. Instead of parsing them, this script reuses tssdistance.py's
nearest-TSS-by-genomic-coordinate approach and extends it to also return
the matched gene's identity (locus_tag/gene/coordinates/strand), not just
the distance.
"""
import os
import sys

import numpy as np
import pandas as pd

XLSX_PATH = "/workspaces/intern-claude-training/data/reference/ncomms5910_supplementary_data1.nomerge.xlsx"
GFF_PATH = "/workspaces/intern-claude-training/data/reference/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff"
OUT_PATH = "/workspaces/intern-claude-training/data/module5/site_to_gene.tsv"

GFF_COLUMNS = [
    "seqid", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes",
]


def load_binding_sites(xlsx_path=XLSX_PATH):
    """Load the nomerge binding-site table and coerce '# of motifs' to numeric.

    header=1 skips the paper's title row above the real header. Rows where
    '# of motifs' is '-' (unknown) are dropped -- they can't be used for
    either hypothesis, which is both keyed on this column.
    """
    if not os.path.isfile(xlsx_path):
        sys.exit(f"Error: binding-site table not found: {xlsx_path}")
    df = pd.read_excel(xlsx_path, header=1)
    df = df.dropna(subset=["Peak"]).reset_index(drop=True)
    df["n_motifs"] = pd.to_numeric(df["# of motifs"], errors="coerce")
    n_dropped = df["n_motifs"].isna().sum()
    df = df.dropna(subset=["n_motifs"]).reset_index(drop=True)
    df["n_motifs"] = df["n_motifs"].astype(int)
    print(f"Binding sites loaded: {len(df)} usable (dropped {n_dropped} with unknown '# of motifs')")
    return df


def load_genes(gff_path=GFF_PATH):
    """Load the annotation GFF and extract locus_tag/gene from the attributes column."""
    if not os.path.isfile(gff_path):
        sys.exit(f"Error: annotation GFF not found: {gff_path}")
    df = pd.read_csv(gff_path, sep="\t", header=None, names=GFF_COLUMNS)
    df["locus_tag"] = df["attributes"].str.extract(r"locus_tag=([^;]+)")
    df["gene"] = df["attributes"].str.extract(r"gene=([^;]+)")
    df["tss"] = np.where(df["strand"] == "+", df["start"], df["end"]).astype(float)
    print(f"Genes loaded: {len(df)}")
    return df


def site_midpoints(df):
    """ChIP-exo Start/End midpoint, unrounded -- matches the paper's own .5-valued distances."""
    return (df["ChIP-exo Start"] + df["ChIP-exo End"]) / 2.0


def nearest_gene_indices(site_positions, tss_array):
    """For each site, the index of the nearest gene by TSS distance (broadcasted)."""
    site_positions = np.asarray(site_positions, dtype=float)
    diff = np.abs(site_positions[:, None] - tss_array[None, :])
    return diff.argmin(axis=1), diff.min(axis=1)


def main():
    sites = load_binding_sites()
    regulatory = sites[sites["Location"] == "regulatory"].reset_index(drop=True)
    print(f"Regulatory-location sites (candidates for gene mapping): {len(regulatory)}")

    genes = load_genes()
    midpoints = site_midpoints(regulatory)
    nearest_idx, nearest_dist = nearest_gene_indices(midpoints, genes["tss"].to_numpy())

    matched = genes.iloc[nearest_idx].reset_index(drop=True)
    out = pd.DataFrame({
        "Peak": regulatory["Peak"].to_numpy(),
        "Transcription_Unit": regulatory["Transcription Unit"].to_numpy(),
        "n_motifs": regulatory["n_motifs"].to_numpy(),
        "S_N_ratio": regulatory["S/N ratio"].to_numpy(),
        "paper_distance_to_tss": regulatory["Distance to TSS"].to_numpy(),
        "matched_locus_tag": matched["locus_tag"].to_numpy(),
        "matched_gene": matched["gene"].to_numpy(),
        "matched_start": matched["start"].to_numpy(),
        "matched_end": matched["end"].to_numpy(),
        "matched_strand": matched["strand"].to_numpy(),
        "computed_distance_bp": nearest_dist,
    })

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"Site-to-gene mapping written to: {OUT_PATH} ({len(out)} rows)")

    dup_genes = out["matched_gene"][out["matched_gene"].duplicated()].unique()
    if len(dup_genes) > 0:
        print(f"Note: {len(dup_genes)} genes are targeted by more than one peak (expected, e.g. multi-site operons): {sorted(dup_genes)}")


if __name__ == "__main__":
    main()
