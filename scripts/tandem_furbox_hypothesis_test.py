#!/usr/bin/env python3
"""Tandem Fur-box 가설에 대해 H1(n vs S/N ratio)과 H2(n vs log2FC)를 검정한다.

H1: tandem Fur-box 개수 n이 늘어날수록 ChIP-exo S/N ratio가 증가한다 (모든
    usable binding site 대상 -- occupancy/signal은 site 자체의 물리적 성질이며
    genomic location과는 무관하다는 가정).
H2: n이 늘어나도 해당 gene의 iron-replete 대비 iron-starved log2FC는
    비례해서 커지지 않는다 -- holo-Fur dimer 하나만으로도 이미 RNAP/sigma를
    완전히 가릴 수 있으므로 saturate될 것으로 예상.

두 Spearman correlation(H1의 rho vs H2의 rho)을 비교하는 것이 핵심 결과다:
n<->log2FC 관계가 n<->S/N ratio 관계보다 뚜렷하게 약하거나 평평하면 H2가
지지된다.
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
    groups = [df.loc[df[n_col] == n, value_col].dropna() for n in sorted(df[n_col].unique())]   # n=1,2,3,4인 행들의 value_col 값들을 그룹화하여 리스트로 모은다. dropna()로 결측치 제거.
    labels = sorted(df[n_col].unique()) # 그림 x축에 "1, 2, 3, 4"라고 표시할 눈금 이름표
    fig, ax = plt.subplots(figsize=(6, 5))  # 그래프 크기 설정
    ax.boxplot(groups, tick_labels=labels)  # boxplot 그리기
    for i, g in enumerate(groups, start=1): # enumerate()로 그룹 인덱스(i)와 그룹 값(g)을 동시에 순회.
        jitter = 0.08 * (pd.Series(range(len(g))) % 2 - 0.5).to_numpy() # jitter를 만들어서 점들을 좌우로 살짝 흩뿌려서 겹치지 않게 한다. (0.08은 흩뿌리는 정도)
        ax.scatter([i] * len(g) + jitter, g, alpha=0.5, color="tab:blue", s=15) # 각 그룹의 점들을 scatter로 찍는다. x좌표는 그룹 인덱스(i)에 jitter를 더해서 흩뿌리고, y좌표는 그룹 값(g)이다. alpha=0.5로 반투명하게 하고, color="tab:blue"로 파란색으로 표시한다.
    ax.set_xlabel("n (# of tandem Fur boxes)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Figure written to: {out_path}")


def run_correlation(df, n_col, value_col, label, unit="sites"):
    sub = df.dropna(subset=[n_col, value_col])  # n_col이나 value_col 값이 결측치 (NaN)인 행(상관관계 계산 불가)을 제거.(defensive coding)
    rho, p = stats.spearmanr(sub[n_col], sub[value_col])    # Spearman 상관계수 rho와 그 유의성(p-value)을 한 번에 계산. rho: n과 value(S/N ratio 또는 log2FC) 간의 monotonic(단조) 관계를 평가.
    group_vals = [sub.loc[sub[n_col] == n, value_col] for n in sorted(sub[n_col].unique())] # n_col 값이 같은 행들을 그룹화하여 value_col 값들을 리스트로 모은다. 예: [n=1인 행들의 value_col 값들, n=2인 행들의 value_col 값들, ...]
    h_stat, kruskal_p = stats.kruskal(*group_vals)  # *groups_vals로 리스트를 풀어서 stats.kruskal(group1, group2, group3, group4) 형태로 넘긴다. Kruskal-Wallis 함수는 그룹 간의 중앙값 차이를 평가하는 비모수적 검정이다. H0: 모든 그룹의 중앙값이 같다.
    print(f"\n--- {label} ---")
    print(f"n = {len(sub)} {unit}, groups (n=1..4) sizes: {[len(g) for g in group_vals]}")  # 전체 표본 크기와 단위를 출력하고, 각 그룹(n=1,2,3,4)의 크기를 리스틀로 출력한다.
    print(f"Spearman rho = {rho:.3f}, p = {p:.4f}") # rho와 p-value를 소수점 3자리와 4자리로 출력한다.
    print(f"Kruskal-Wallis H = {h_stat:.3f}, p = {kruskal_p:.4f}")  # Kruskal-Wallis 검정 통계량과 p-value를 소수점 3자리와 4자리로 출력한다.
    print("Group medians:")
    for n, g in zip(sorted(sub[n_col].unique()), group_vals):   # 각 그룹의 중앙값, 평균, 표본 크기를 출력한다. zip()을 사용하여 n 값과 해당 그룹의 value_col 값들을 동시에 순회한다.
        print(f"  n={n}: median={g.median():.3f}, mean={g.mean():.3f} (N={len(g)})")
    return rho, p   # rho와 p만 반환


def main():
    # H1: 모든 usable binding site 대상, n vs S/N ratio
    sites = load_binding_sites()
    h1_rho, h1_p = run_correlation(sites, "n_motifs", "S/N ratio", "H1: n vs S/N ratio (all sites)")
    boxplot_by_n(
        sites, "n_motifs", "S/N ratio",
        ylabel="ChIP-exo S/N ratio",
        title=f"H1: n vs S/N ratio (rho={h1_rho:.2f}, p={h1_p:.3g})",
        out_path=FIG1_PATH,
    )

    # H2: gene에 매핑된 regulatory-location site 대상, n vs log2FC
    if not os.path.isfile(HYPOTHESIS_TABLE_PATH):
        print(
            f"\n{HYPOTHESIS_TABLE_PATH} not found -- run scripts/tandem_furbox_expression.py "
            "first (needs SRR1168134/SRR1168136 downloaded+aligned and all 4 *_linear.gff built)."
        )
        return
    combined = pd.read_csv(HYPOTHESIS_TABLE_PATH, sep="\t")
    # combined는 peak(site) 단위라, gene 하나에 peak가 여러 개면 그 gene의 log2FC가
    # 행마다 그대로 복제되어 있다 (pseudoreplication). 상관관계 검정 전에 gene당
    # 한 행만 남기도록 집계한다: n_motifs는 그 gene의 site들 중 가장 큰 tandem
    # Fur-box 개수(가장 강한 결합 클러스터)를 대표값으로 쓰고, log2FC는 gene당
    # 이미 유일한 값이므로 첫 값을 그대로 쓴다.
    combined_by_gene = combined.groupby("matched_locus_tag", as_index=False).agg(   # mathced_locus_tag 기준으로 그룹화하여 gene 단위로 집계한다. as_index=False로 그룹화된 열을 인덱스로 따로 만들지 않는다.
        matched_gene=("matched_gene", "first"), # matched_gene 열의 첫 번째 값을 대표값으로 사용한다. (gene 이름)
        n_motifs=("n_motifs", "max"),   # n_motifs 열의 최대값을 대표값으로 사용한다. (gene에 매핑된 site들 중 가장 큰 tandem Fur-box 개수)
        log2FC=("log2FC", "first"), # log2FC 열의 첫 번째 값을 대표값으로 사용한다. (gene의 iron-starved 대비 iron-replete log2FC)
    )
    h2_rho, h2_p = run_correlation(
        combined_by_gene, "n_motifs", "log2FC", "H2: n vs log2FC (one row per gene)", unit="genes"  # unit="genes"로 gene 단위로 상관관계 검정한다.
    )
    boxplot_by_n(   # gene 단위로 boxplot 그리기
        combined_by_gene, "n_motifs", "log2FC",
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
