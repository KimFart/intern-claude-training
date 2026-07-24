#!/usr/bin/env python3
"""Fur ChIP-exo 결합 부위 표를 GFF 트랙으로 변환한다 (Module 4, Exercise 8).

MetaScope는 엑셀 표를 읽지 못하므로, iron-replete 결합 부위(143개)를
data/reference/fur_sites.gff로 저장해 다른 트랙(annotation, rnaseq.gff)과
나란히 시각화할 수 있게 한다.

컬럼 1은 다른 트랙과 맞추기 위해 버전 접미사 없는 'NC_000913'을 그대로 쓰고,
score 컬럼에는 S/N ratio를 넣어 결합 강도를 시각적으로 구분할 수 있게 한다.
"""
import os
import sys

import pandas as pd

XLSX_PATH = "/workspaces/intern-claude-training/data/reference/ncomms5910_supplementary_data1.xlsx"
OUT_GFF_PATH = "/workspaces/intern-claude-training/data/reference/fur_sites.gff"


def load_binding_sites(xlsx_path=XLSX_PATH):
    """iron-replete Fur 결합 부위 표를 읽어 DataFrame으로 반환한다.

    엑셀 파일은 실제 헤더 위에 제목 행이 있으므로 header=1로 읽는다. 한 행은 전부
    NaN인 placeholder이므로 dropna(subset=["Peak"])로 제거한다. Binding Conditiona
    열에 'R'이 포함된 행('R' 또는 'R/S')이 iron-replete 부위이며, 이는 Exercise 5/7과
    동일한 필터링 기준이다.
    """
    if not os.path.isfile(xlsx_path):
        sys.exit(f"Error: binding-site table not found: {xlsx_path}")
    df = pd.read_excel(xlsx_path, header=1)
    df = df.dropna(subset=["Peak"])
    df = df[df["Binding Conditiona"].astype(str).str.contains("R")]
    print(f"Iron-replete binding sites: {len(df)}")
    return df.reset_index(drop=True)


def build_gff_line(row):
    """결합 부위 한 행을 GFF 한 줄로 변환한다.

    Transcription Unit이 없는 행(NaN)도 그대로 남긴다 -- Exercise 8은 iron-replete
    부위 143개 전부가 fur_sites.gff에 한 줄씩 있는지 확인하라고 했으므로, 유전자
    이름이 없다고 행을 빼면 안 된다. str(NaN)은 그대로 'nan' 문자열이 되는데, 이는
    sitesformeme.py의 positive_control()에서 이미 쓰고 있는 방식과 같다.
    """
    start = int(row["ChIP-exo Start"])
    end = int(row["ChIP-exo End"])
    score = row["S/N ratio"]
    name = str(row["Transcription Unit"])
    return f"NC_000913\tfurchip\tfur_site\t{start}\t{end}\t{score}\t.\t.\tname={name}\n"


def main():
    df = load_binding_sites()
    with open(OUT_GFF_PATH, "w") as gff:
        for _, row in df.iterrows():
            gff.write(build_gff_line(row))
    print(f"GFF rows written: {len(df)}")
    print(f"Output written to: {OUT_GFF_PATH}")


if __name__ == "__main__":
    main()
