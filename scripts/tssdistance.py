#!/usr/bin/env python3
"""Fur binding site(결합 부위)와 가장 가까운 TSS 사이의 거리를 계산한다.

두 개의 데이터 소스를 사용한다:
- 결합 부위 좌표: ncomms5910_supplementary_data1.xlsx (Seo et al. 2014 Supplementary Data 1,
  Exercise 5와 동일하게 iron-replete 부위만 사용)
- annotated gene: ec_annotation_20100903_DHK_cSRNA_with_ortho.gff (EcoCyc 14.1, 모든 행이
  'gene' feature이므로 strand에 따라 start/end에서 TSS를 직접 유도한다)

"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

XLSX_PATH = "/workspaces/intern-claude-training/data/reference/ncomms5910_supplementary_data1.xlsx"
GFF_PATH = "/workspaces/intern-claude-training/data/reference/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff"
HIST_PATH = "/workspaces/intern-claude-training/data/reference/fur_tss_distance_histogram.png"
TSS_DISTANCE_THRESHOLD = 500  # bp

GFF_COLUMNS = [
    "seqid", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes",
]


def load_binding_sites(xlsx_path=XLSX_PATH):
    """iron-replete Fur 결합 부위 표를 읽어 DataFrame으로 반환한다.

    엑셀 파일은 실제 헤더 위에 제목 행이 있으므로 header=1로 읽는다. 한 행은 전부
    NaN인 placeholder이므로 dropna(subset=["Peak"])로 제거한다. Binding Condition
    열에 'R'이 포함된 행('R' 또는 'R/S')이 iron-replete 부위이며, 이는 Exercise 5
    와 동일한 필터링 기준이다.
    """
    if not os.path.isfile(xlsx_path):
        sys.exit(f"Error: binding-site table not found: {xlsx_path}")
    df = pd.read_excel(xlsx_path, header=1)
    df = df.dropna(subset=["Peak"])
    df = df[df["Binding Conditiona"].astype(str).str.contains("R")]
    print(f"Iron-replete binding sites: {len(df)}")
    return df.reset_index(drop=True)


def load_gene_tss(gff_path=GFF_PATH):
    """GFF annotation을 읽어 각 유전자의 TSS 좌표 배열을 반환한다.

    이 파일의 모든 행의 feature 컬럼은 'gene'이므로, 각 행에서 strand가 '+'이면 
    start를, '-'이면 end를 TSS로 사용한다. 이 파일은 고정된 입력이므로 다른 
    feature 타입이 섞여 있는지 여부는 print로만 확인하고 (경고), 실행을 막지는 않는다.
    """
    if not os.path.isfile(gff_path):
        sys.exit(f"Error: annotation GFF not found: {gff_path}")
    df = pd.read_csv(gff_path, sep="\t", header=None, names=GFF_COLUMNS)
    non_gene = df[df["feature"] != "gene"]
    if len(non_gene) > 0:
        print(
            f"Note: {len(non_gene)} rows are not 'gene' features and were included "
            "anyway using the same start/end-by-strand rule."
        )
    tss = np.where(df["strand"] == "+", df["start"], df["end"]).astype(float)
    print(f"Genes loaded: {len(df)}")
    return tss


def site_midpoints(df):
    """ChIP-exo Start/End의 중간점을 반올림 없이 그대로 계산한다.

    논문의 Distance to TSS 열에 47.5, 11.5 같은 .5 값이 있는 것으로 보아, 원 저자들도
    정수로 반올림하지 않은 (start+end)/2를 그대로 사용했음을 알 수 있다. 따라서 여기서도
    동일한 방식(반올림 없는 float 평균)을 사용해야 cross-check가 의미를 가진다.
    """
    return (df["ChIP-exo Start"] + df["ChIP-exo End"]) / 2.0


def nearest_tss_distance(site_positions, tss_array):
    """각 결합 부위(site)에서 가장 가까운 전사 시작 위치(TSS)까지의 거리(bp)를 계산한다.

    [동작 원리]
    1. site_positions[:, None]  -> (N, 1) 세로 배열로 변경
    2. tss_array[None, :]       -> (1, M) 가로 배열로 변경
    3. 브로드캐스팅 뺄셈            -> (N, M) 모든 조합의 거리 행렬 생성
    4. min(axis=1)              -> 각 site(행) 기준 가장 가까운 거리 선택

    ※ N=143, M=4595 기준 약 65만 번 연산으로, NumPy 브로드캐스팅 시 1ms 미만 소요.
    """    
    site_positions = np.asarray(site_positions, dtype=float)
    diff = np.abs(site_positions[:, None] - tss_array[None, :])
    return diff.min(axis=1)


def cross_check(computed_distances, df):
    """스크립트가 계산한 TSS 거리와 논문(Supplementary Data 1)의 기존 거리를 비교 검증한다.
    
    [검증 방식 및 예외 처리]
    1. 데이터 정제: 논문 데이터의 결측치/문자열('-')을 NaN으로 처리 후 숫자 행만 추출
    2. 절댓값 비교: 방향(부호, +/-) 정의가 논문과 다를 수 있어 순수 거리만 비교
    3. 허용 가능한 경고: '기계적 최단 거리' vs '생물학적 특정 유전자 거리'의
    정의 차이로 오차가 생길 수 있으므로, 오차가 커도 프로그램 중단 없이 경고만 출력
    """
    # 1. 'Distance to TSS' 열의 문자열('-' 등)을 NaN으로 바꾸고 숫자 행만 골라냄
    paper_distance = pd.to_numeric(df["Distance to TSS"], errors="coerce")
    mask = paper_distance.notna()
    
    # 2. 비교할 유효 데이터가 없으면 안내 후 종료
    if mask.sum() == 0:
        print("Cross-check: no numeric 'Distance to TSS' values found in the table, skipping.")
        return
    
    # 3. 부호 기준 차이를 방지하기 위해 절댓값끼리 차이(bp) 계산
    diff = np.abs(computed_distances[mask.to_numpy()]) - paper_distance[mask].abs().to_numpy()
    mean_abs_offset = np.abs(diff).mean()
    print(
        f"Cross-check against paper's 'Distance to TSS' column "
        f"({mask.sum()} sites with a numeric value):"
    )
    print(f"  mean |computed - paper| offset: {mean_abs_offset:.1f} bp")
    if mean_abs_offset > 100:
        print(
            "  Note: the average offset is fairly large. This can be expected if the "
            "paper's distance is to a specific assigned Transcription Unit rather than "
            "the literal nearest gene genome-wide -- worth a quick look, not necessarily a bug."
        )


def plot_histogram(distances, out_path=HIST_PATH):
    """결합 부위 - 가장 가까운 TSS 거리 분포를 히스토그램으로 그려 저장한다."""
    plt.figure(figsize=(8, 5))
    plt.hist(distances, bins=30, edgecolor="black")
    plt.xlabel("Distance to nearest TSS (bp)")
    plt.ylabel("Number of Fur binding sites")
    plt.title("Fur binding sites: distance to nearest TSS")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"Histogram written to: {out_path}")


def main():
    sites_df = load_binding_sites()
    tss_array = load_gene_tss()
    midpoints = site_midpoints(sites_df)
    distances = nearest_tss_distance(midpoints, tss_array)

    print(f"Minimum distance to nearest TSS: {distances.min():.1f} bp")
    print(f"Median distance to nearest TSS: {np.median(distances):.1f} bp")
    within_threshold = (distances <= TSS_DISTANCE_THRESHOLD).sum()
    fraction = within_threshold / len(distances)
    print(
        f"Fraction within {TSS_DISTANCE_THRESHOLD} bp of a TSS: "
        f"{fraction:.3f} ({within_threshold}/{len(distances)})"
    )

    cross_check(distances, sites_df)
    plot_histogram(distances)


if __name__ == "__main__":
    main()
