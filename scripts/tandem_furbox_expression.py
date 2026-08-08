#!/usr/bin/env python3
"""Tandem Fur-box 가설(H2)을 위한 gene 단위 log2FC 테이블을 만든다.

site_to_gene.py의 매핑 결과에 있는 각 Fur-regulated gene에 대해, 4개의
*_linear.gff 트랙(notebooks/makegff.py를 --log_scale 없이 실행해서 만든 것 --
depth를 선형(linear)으로 그대로 합산할 수 있어야 하므로, 기존 --log_scale GFF는
이 용도로 쓸 수 없다는 점은 plan의 note 참고)에서 gene body 구간의 per-base
RNA-seq depth를 합산한다. 그 다음 각 replicate를 자신의 total mapped-read count로
정규화하고, condition별로 두 replicate를 평균낸 뒤 log2FC(starve / replete)를 계산한다.

*_linear.gff의 각 행은 depth가 0이 아닌 위치 하나씩을 나타낸다(start == end,
염기 하나당 한 행 -- notebooks/makegff.py는 depth가 여러 염기에 걸쳐 똑같이 
유지되는 구간을 한 행으로 합치는 run-length encoding을 하지 않음). 따라서
gene body의 총 depth는, gene과 strand가 같고 위치가 그 gene의 [start, end]
범위 안에 들어오는 행들만 골라서 그 score 값을 다 더하면 구한다.
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

DATA_MOD5 = "/workspaces/intern-claude-training/data/module5"
SITE_TO_GENE_PATH = os.path.join(DATA_MOD5, "site_to_gene.tsv")
OUT_PATH = os.path.join(DATA_MOD5, "tandem_furbox_hypothesis_table.tsv")

# SRR -> (condition, replicate 번호)
REPLICATES = {
    "SRR1168133": ("replete", 1),
    "SRR1168134": ("replete", 2),
    "SRR1168135": ("starve", 1),
    "SRR1168136": ("starve", 2),
}
PSEUDOCOUNT = 1.0
CPM_SCALE = 1e6


def total_mapped_reads(srr):
    """sorted BAM의 total mapped reads. samtools idxstats의 3번째 열 합계로 구한다."""
    bam_path = os.path.join(DATA_MOD5, f"{srr}_sorted.bam")
    if not os.path.isfile(bam_path):
        sys.exit(f"Error: sorted BAM not found: {bam_path}")
    result = subprocess.run(
        ["samtools", "idxstats", bam_path], capture_output=True, text=True, check=True
    )
    total = sum(int(line.split("\t")[2]) for line in result.stdout.strip().splitlines())
    return total


def load_linear_gff(srr):
    """*_linear.gff를 읽는다 (per-base depth, depth가 0이 아닌 위치마다 한 행)."""
    gff_path = os.path.join(DATA_MOD5, f"{srr}_linear.gff")
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
    """gene별 depth 합계 (strand가 일치하고 위치가 [start, end] 안에 있는 것만)."""
    sums = {}
    for strand in ("+", "-"):                                   # strand별로 따로 처리
        strand_rows = gff_df[gff_df["strand"] == strand]        # strand가 일치하는 행만 남김
        positions = strand_rows["start"].to_numpy()
        scores = strand_rows["score"].to_numpy()
        for _, gene in genes[genes["matched_strand"] == strand].iterrows():
            mask = (positions >= gene["matched_start"]) & (positions <= gene["matched_end"])    # boolean mask
            sums[gene["matched_locus_tag"]] = float(scores[mask].sum())                         # 더하기만 함
    return sums


def main():
    if not os.path.isfile(SITE_TO_GENE_PATH):
        sys.exit(f"Error: {SITE_TO_GENE_PATH} not found. Run scripts/site_to_gene.py first.")
    sites = pd.read_csv(SITE_TO_GENE_PATH, sep="\t")

    genes = sites.drop_duplicates(subset=["matched_locus_tag"])[                                    # matched_locus_tag 기준으로 중복 제거
        ["matched_locus_tag", "matched_gene", "matched_start", "matched_end", "matched_strand"]     # 제거한 DataFrame에서 gene 자체에 대한 정보(locus_tag, gene 이름, 좌표, strand)만 남긴다.
    ]
    print(f"Unique genes to quantify: {len(genes)} (from {len(sites)} peaks)")

    normalized = {}                                                                                 # locus_tag -> {"replete": [rep1, rep2], "starve": [rep1, rep2]} 형태로 채워짐
    for srr, (condition, rep) in REPLICATES.items():
        print(f"Processing {srr} ({condition} rep{rep})...")
        mapped_reads = total_mapped_reads(srr)                                                      # 각 SRR에 대해 total mapped reads를 구한다.
        gff_df = load_linear_gff(srr)                                                               # 각 SRR에 대한 linear GFF 파일을 로드한다.
        depth_sums = gene_depth_sums(gff_df, genes)                                                 # 각 gene에 대해 depth 합계를 계산한다.
        for locus_tag, depth_sum in depth_sums.items():                                             # 각 gene에 대해 CPM을 계산하고 normalized 딕셔너리에 저장한다.
            cpm = depth_sum * CPM_SCALE / mapped_reads                                              # CPM 계산: (gene depth 합계 / total mapped reads) * 1e6
            normalized.setdefault(locus_tag, {"replete": [], "starve": []})[condition].append(cpm)  # condition에 따라 replete 또는 starve 리스트에 CPM 값 추가
        del gff_df                                                                                  # 메모리 절약을 위해 gff_df 삭제

    rows = []   # gene 단위 log2FC 계산 결과를 담을 리스트. 최종적으로 gene 하나당 dict 하나씩 61개가 쌓이고, DataFrame으로 변환된다.
    for _, gene in genes.iterrows():    # gene 단위로 반복하면서 log2FC 계산
        locus_tag = gene["matched_locus_tag"]   # locus_tag를 기준으로 normalized 딕셔너리에서 replete와 starve replicate CPM 값을 가져온다.
        replete_reps = normalized[locus_tag]["replete"] # replete 조건의 replicate CPM 값들을 가져온다.
        starve_reps = normalized[locus_tag]["starve"]   # starve 조건의 replicate CPM 값들을 가져온다.
        replete_mean = float(np.mean(replete_reps)) # replete 조건의 replicate CPM 값들의 평균을 계산한다.
        starve_mean = float(np.mean(starve_reps))   # starve 조건의 replicate CPM 값들의 평균을 계산한다.
        log2fc = np.log2((starve_mean + PSEUDOCOUNT) / (replete_mean + PSEUDOCOUNT)) # log2FC 계산: log2((starve_mean + pseudocount) / (replete_mean + pseudocount)), pseudocount를 더해 0으로 나누거나 log0이 되는 것을 방지한다.
        rows.append({   # gene 단위 log2FC 계산 결과를 담은 dict를 rows 리스트에 추가한다.
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
    gene_table = pd.DataFrame(rows) # 61개 gene에 대한 dict 리스트를 하나의 DataFrame으로 변환한다.

    combined = sites.merge(gene_table, on=["matched_locus_tag", "matched_gene"], how="left")    # site 단위 정보(sites)와 gene 단위 log2FC 계산 결과(gene_table)를 locus_tag와 gene 이름을 기준으로 병합한다. how="left"로 하여 sites에 있는 모든 peak 정보를 유지하고, gene_table에 없는 경우 NaN으로 채운다.
    os.makedirs(DATA_MOD5, exist_ok=True)
    combined.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"Hypothesis table written to: {OUT_PATH} ({len(combined)} rows, {len(gene_table)} unique genes)")


if __name__ == "__main__":
    main()
