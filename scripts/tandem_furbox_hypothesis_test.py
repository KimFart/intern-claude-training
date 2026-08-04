#!/usr/bin/env python3
"""Test H1 (n vs S/N ratio) and H2 (n vs log2FC) for the tandem Fur-box hypothesis.

H1: as tandem Fur-box count n increases, ChIP-exo S/N ratio increases (all
    usable binding sites -- occupancy/signal is a physical property of the
    site itself, independent of genomic location).
H2: as n increases, the regulated gene's iron-replete-vs-starved log2FC does
    NOT scale proportionally -- expected to saturate, since one holo-Fur
    dimer can already fully occlude RNAP/sigma.

Comparing the two Spearman correlations (H1's rho vs H2's rho) is the core
result: H2 is supported if the n<->log2FC relationship is markedly weaker/
flatter than the n<->S/N ratio relationship.
"""
import os

import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

from site_to_gene import load_binding_sites

DATA_MOD5 = "/workspaces/intern-claude-training/data/module5"
HYPOTHESIS_TABLE_PATH = os.path.join(DATA_MOD5, "tandem_furbox_hypothesis_table.tsv")
FIG1_PATH = os.path.join(DATA_MOD5, "fig1_n_vs_sn_ratio.png")
FIG2_PATH = os.path.join(DATA_MOD5, "fig2_n_vs_log2fc.png")


def boxplot_by_n(df, n_col, value_col, ylabel, title, out_path):
    groups = [df.loc[df[n_col] == n, value_col].dropna() for n in sorted(df[n_col].unique())]
    labels = sorted(df[n_col].unique())
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.boxplot(groups, tick_labels=labels)
    for i, g in enumerate(groups, start=1):
        jitter = 0.08 * (pd.Series(range(len(g))) % 2 - 0.5).to_numpy()
        ax.scatter([i] * len(g) + jitter, g, alpha=0.5, color="tab:blue", s=15)
    ax.set_xlabel("n (# of tandem Fur boxes)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Figure written to: {out_path}")


def run_correlation(df, n_col, value_col, label):
    sub = df.dropna(subset=[n_col, value_col])
    rho, p = stats.spearmanr(sub[n_col], sub[value_col])
    group_vals = [sub.loc[sub[n_col] == n, value_col] for n in sorted(sub[n_col].unique())]
    h_stat, kruskal_p = stats.kruskal(*group_vals)
    print(f"\n--- {label} ---")
    print(f"n = {len(sub)} sites, groups (n=1..4) sizes: {[len(g) for g in group_vals]}")
    print(f"Spearman rho = {rho:.3f}, p = {p:.4f}")
    print(f"Kruskal-Wallis H = {h_stat:.3f}, p = {kruskal_p:.4f}")
    print("Group medians:")
    for n, g in zip(sorted(sub[n_col].unique()), group_vals):
        print(f"  n={n}: median={g.median():.3f}, mean={g.mean():.3f} (N={len(g)})")
    return rho, p


def main():
    # H1: all usable binding sites, n vs S/N ratio
    sites = load_binding_sites()
    h1_rho, h1_p = run_correlation(sites, "n_motifs", "S/N ratio", "H1: n vs S/N ratio (all sites)")
    boxplot_by_n(
        sites, "n_motifs", "S/N ratio",
        ylabel="ChIP-exo S/N ratio",
        title=f"H1: n vs S/N ratio (rho={h1_rho:.2f}, p={h1_p:.3g})",
        out_path=FIG1_PATH,
    )

    # H2: regulatory-location sites mapped to a gene, n vs log2FC
    if not os.path.isfile(HYPOTHESIS_TABLE_PATH):
        print(
            f"\n{HYPOTHESIS_TABLE_PATH} not found -- run scripts/tandem_furbox_expression.py "
            "first (needs SRR1168134/SRR1168136 downloaded+aligned and all 4 *_linear.gff built)."
        )
        return
    combined = pd.read_csv(HYPOTHESIS_TABLE_PATH, sep="\t")
    h2_rho, h2_p = run_correlation(combined, "n_motifs", "log2FC", "H2: n vs log2FC (regulatory sites mapped to a gene)")
    boxplot_by_n(
        combined, "n_motifs", "log2FC",
        ylabel="log2FC (iron-starved / iron-replete)",
        title=f"H2: n vs log2FC (rho={h2_rho:.2f}, p={h2_p:.3g})",
        out_path=FIG2_PATH,
    )

    print("\n--- H1 vs H2 comparison ---")
    print(f"H1 (n vs S/N ratio):  rho = {h1_rho:.3f}, p = {h1_p:.4f}")
    print(f"H2 (n vs log2FC):     rho = {h2_rho:.3f}, p = {h2_p:.4f}")
    if abs(h2_rho) < abs(h1_rho):
        print("H2's correlation is weaker than H1's, consistent with a saturating repression response.")
    else:
        print("H2's correlation is NOT weaker than H1's -- saturation hypothesis is not supported by this comparison.")


if __name__ == "__main__":
    main()
