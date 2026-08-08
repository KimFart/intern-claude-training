#!/usr/bin/env python3
"""Fur 결합 site를 해당 전사단위(TU/operon)의 유전자에 안전하게 매핑한다.

기존의 단순 nearest-TSS 방식은 peak 바로 옆에 있는 '다른 operon'의 유전자를
잘못 선택할 수 있다. 이 스크립트는 논문 표의 Transcription Unit(TU)을 먼저
annotation의 gene symbol들로 분해하여, TU에 속하는 유전자만 후보로 제한한다.
따라서 인접하지만 다른 operon에 속한 유전자가 더 가까워도 선택하지 않는다.

중요한 설계 원칙
----------------
1) TU를 해석할 수 없거나 해석이 여러 개로 모호하면, genome 전체의 nearest gene으로
   억지 매핑하지 않는다. 대신 mapping_status에 이유를 기록하고 matched_* 값은 비운다.
2) `Distance to TSS`가 논문 표에 있으면, TU 내부 후보 TSS 중 그 거리와 가장 잘 맞는
   TSS를 선택한다. 이는 internal promoter/TSS가 있는 TU도 일부 반영한다.
3) 예외적인 TU는 `tu_primary_gene_map.tsv`의 수동 curated override로 해결한다.
   이 파일은 권장 입력이며, 최소 열은 `transcription_unit`와 `locus_tag`이다.

필수 입력
---------
- ncomms5910_supplementary_data1.nomerge.xlsx
- ec_annotation_20100903_DHK_cSRNA_with_ortho.gff

선택 입력
---------
- data/module5/tu_primary_gene_map.tsv
  예시:
  transcription_unit\tlocus_tag
  fhuACDB\tb0150
  entCEBAH\tb0585

출력
----
- data/module5/site_to_gene.tsv
"""

import os
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

DATA_REF = "/workspaces/intern-claude-training/data/reference"
DATA_MOD5 = "/workspaces/intern-claude-training/data/module5"
XLSX_PATH = os.path.join(DATA_REF, "ncomms5910_supplementary_data1.nomerge.xlsx")
GFF_PATH = os.path.join(DATA_REF, "ec_annotation_20100903_DHK_cSRNA_with_ortho.gff")
OVERRIDE_PATH = os.path.join(DATA_MOD5, "tu_primary_gene_map.tsv")
OUT_PATH = os.path.join(DATA_MOD5, "site_to_gene.tsv")

GFF_COLUMNS = [
    "seqid", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes",
]


def normalize_tu(value):
    """TU/gene 이름 비교용 표준화: 대소문자 차이와 기호를 제거한다."""
    if pd.isna(value):
        return ""
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def load_binding_sites(xlsx_path=XLSX_PATH):
    """논문 supplementary Excel을 읽고 motif 수가 알려진 peak만 남긴다."""
    if not os.path.isfile(xlsx_path):
        sys.exit(f"Error: binding-site table not found: {xlsx_path}")

    # Excel 첫 행은 표 제목이고, 두 번째 행이 실제 column header인 형식이다.
    sites = pd.read_excel(xlsx_path, header=1)
    sites = sites.dropna(subset=["Peak"]).reset_index(drop=True)
    sites["n_motifs"] = pd.to_numeric(sites["# of motifs"], errors="coerce")
    n_unknown = sites["n_motifs"].isna().sum()
    sites = sites.dropna(subset=["n_motifs"]).reset_index(drop=True)
    sites["n_motifs"] = sites["n_motifs"].astype(int)
    print(f"Binding sites loaded: {len(sites)} usable (dropped {n_unknown} with unknown motif count)")
    return sites


def load_genes(gff_path=GFF_PATH):
    """annotation GFF에서 실제 gene feature만 읽고 TSS를 계산한다."""
    if not os.path.isfile(gff_path):
        sys.exit(f"Error: annotation GFF not found: {gff_path}")

    genes = pd.read_csv(gff_path, sep="\t", header=None, names=GFF_COLUMNS, comment="#")

    # CDS, mRNA, region 등 다른 feature가 nearest-TSS 후보가 되는 것을 막는다.
    genes = genes[genes["feature"] == "gene"].copy()
    genes["locus_tag"] = genes["attributes"].str.extract(r"locus_tag=([^;]+)")
    genes["gene"] = genes["attributes"].str.extract(r"gene=([^;]+)")
    genes = genes.dropna(subset=["locus_tag", "gene"])
    genes = genes[genes["strand"].isin(["+", "-"])].copy()

    # + strand에서는 작은 좌표(start), - strand에서는 큰 좌표(end)가 5' TSS이다.
    genes["tss"] = np.where(genes["strand"] == "+", genes["start"], genes["end"]).astype(float)
    genes["gene_key"] = genes["gene"].map(normalize_tu)
    genes = genes.reset_index(drop=True)
    print(f"Gene features loaded: {len(genes)}")
    return genes


def load_overrides(override_path=OVERRIDE_PATH):
    """수동 검증한 TU→primary locus_tag override를 읽는다; 없으면 빈 dict를 반환한다."""
    if not os.path.isfile(override_path):
        print(f"Override file not found (optional): {override_path}")
        return {}

    override = pd.read_csv(override_path, sep="\t")
    required = {"transcription_unit", "locus_tag"}
    missing = required - set(override.columns)
    if missing:
        sys.exit(f"Error: override file missing columns: {sorted(missing)}")

    override["tu_key"] = override["transcription_unit"].map(normalize_tu)
    if override["tu_key"].duplicated().any():
        sys.exit("Error: each transcription_unit may appear only once in the override file")
    return dict(zip(override["tu_key"], override["locus_tag"]))


STEM_SUFFIX_PATTERN = re.compile(r"^([a-z]{2,4})([A-Z]{2,})$")


def stem_suffix_segmentation(tu_raw, gene_key_set):
    """'stem + 개별 대문자 접미사' 형태의 operon 표기를 분해한다.

    예: 'fhuACDB' (stem 'fhu' + 접미사 A,C,D,B) → ('fhua', 'fhuc', 'fhud', 'fhub')
    all_tu_segmentations는 TU 문자열이 완전한 gene symbol들의 연결일 때만 동작하므로,
    stem을 공유하고 마지막 글자만 다른 이 표기(예: entCEBAH, fecIR)는 분해하지 못한다.
    재구성한 gene symbol이 annotation에 전부 존재할 때만 결과를 반환하고, 하나라도
    없으면 빈 리스트를 반환해 억지 매핑을 피한다.
    """
    match = STEM_SUFFIX_PATTERN.match(str(tu_raw).strip())
    if not match:
        return []
    stem, suffixes = match.groups()
    gene_symbols = tuple(normalize_tu(stem + letter) for letter in suffixes)
    if any(symbol not in gene_key_set for symbol in gene_symbols):
        return []
    return [gene_symbols]


def all_tu_segmentations(tu_key, gene_keys, max_solutions=20):
    """붙어 있는 TU 문자열을 annotation gene symbol들의 조합으로 모두 분해한다.

    예: 'fepAentD' (하이픈 제거 후) → ('fepa', 'entd')처럼, TU 문자열이 완전한
    gene symbol들을 그대로 이어 붙인 경우를 처리한다. stem+접미사 표기는
    stem_suffix_segmentation이 별도로 처리한다.
    완전탐색이지만 문자열 길이가 짧고 memoization을 사용한다. 분해가 2개 이상이면
    자동 선택은 위험하므로 caller가 ambiguous 상태로 기록한다.
    """
    gene_key_set = set(gene_keys)
    by_first_char = defaultdict(list)
    for key in gene_key_set:
        by_first_char[key[0]].append(key)
    for key_list in by_first_char.values():
        key_list.sort(key=len, reverse=True)

    memo = {}

    def walk(position):
        if position == len(tu_key):
            return [()]
        if position in memo:
            return memo[position]

        solutions = []
        for key in by_first_char.get(tu_key[position], []):
            if not tu_key.startswith(key, position):
                continue
            for tail in walk(position + len(key)):
                solutions.append((key,) + tail)
                if len(solutions) >= max_solutions:
                    memo[position] = solutions
                    return solutions
        memo[position] = solutions
        return solutions

    return walk(0)


def numeric_distance(value):
    """논문 Distance to TSS를 숫자로 바꾼다; 부호가 있더라도 절댓값 비교에 쓴다."""
    try:
        return abs(float(value))
    except (TypeError, ValueError):
        return np.nan


def choose_tu_gene(site, genes, gene_to_rows, gene_keys):
    """한 regulatory peak에 대해 TU 내부의 가장 적절한 TSS gene을 선택한다.

    반환값은 (gene_row 또는 None, status, candidate_gene_symbols, computed_distance).
    논문 거리값이 있으면 peak-TSS 절대거리가 논문 거리와 가장 가까운 후보를 고른다.
    거리값이 없으면 TU의 가장 upstream TSS를 primary gene으로 간주한다.
    """
    tu = site["Transcription Unit"]
    tu_key = normalize_tu(tu)
    if not tu_key:
        return None, "unresolved_missing_TU", "", np.nan

    segmentations = all_tu_segmentations(tu_key, gene_keys)
    if not segmentations:
        # 'fhuACDB'류 stem+접미사 표기는 명시적 gene-name 조합으로는 분해되지 않는다.
        segmentations = stem_suffix_segmentation(tu, gene_to_rows)
    unique_segmentations = {tuple(x) for x in segmentations}
    if len(unique_segmentations) == 0:
        return None, "unresolved_TU_not_parseable", "", np.nan
    if len(unique_segmentations) > 1:
        choices = " | ".join("+".join(x) for x in sorted(unique_segmentations))
        return None, "unresolved_ambiguous_TU_parse", choices, np.nan

    gene_symbols = next(iter(unique_segmentations))
    candidate_rows = []
    for gene_key in gene_symbols:
        candidate_rows.extend(gene_to_rows[gene_key])
    candidates = genes.loc[candidate_rows].copy()
    candidate_text = "+".join(gene_symbols)

    if candidates.empty:
        return None, "unresolved_no_annotated_TU_gene", candidate_text, np.nan

    midpoint = (float(site["ChIP-exo Start"]) + float(site["ChIP-exo End"])) / 2.0
    candidates["peak_tss_distance"] = (candidates["tss"] - midpoint).abs()
    paper_distance = numeric_distance(site.get("Distance to TSS"))

    if not np.isnan(paper_distance):
        # 논문이 제공한 TSS 거리와 가장 잘 맞는 TU 내부 TSS를 선택한다.
        candidates["paper_distance_error"] = (candidates["peak_tss_distance"] - paper_distance).abs()
        selected = candidates.sort_values(["paper_distance_error", "peak_tss_distance", "locus_tag"]).iloc[0]
        status = "mapped_TU_constrained_paper_distance"
    else:
        # 논문 거리 정보가 없으면 operon의 5' 쪽(primary) gene을 사용한다.
        strands = candidates["strand"].unique()
        if len(strands) != 1:
            return None, "unresolved_mixed_strand_TU", candidate_text, np.nan
        if strands[0] == "+":
            selected = candidates.sort_values(["tss", "locus_tag"]).iloc[0]
        else:
            selected = candidates.sort_values(["tss", "locus_tag"], ascending=[False, True]).iloc[0]
        status = "mapped_TU_constrained_primary_TSS"

    return selected, status, candidate_text, float(selected["peak_tss_distance"])


def main():
    sites = load_binding_sites()
    regulatory = sites[sites["Location"].astype(str).str.lower() == "regulatory"].reset_index(drop=True)
    print(f"Regulatory sites to map: {len(regulatory)}")

    genes = load_genes()
    overrides = load_overrides()
    gene_to_rows = defaultdict(list)
    for idx, gene_key in enumerate(genes["gene_key"]):
        gene_to_rows[gene_key].append(idx)
    gene_keys = sorted(gene_to_rows)
    locus_to_row = genes.set_index("locus_tag", drop=False)

    output_rows = []
    for _, site in regulatory.iterrows():
        tu_key = normalize_tu(site["Transcription Unit"])
        selected = None
        candidate_text = ""
        computed_distance = np.nan

        # Curated override가 있으면 가장 신뢰도가 높으므로 자동 규칙보다 우선한다.
        override_locus = overrides.get(tu_key)
        if override_locus is not None:
            if override_locus in locus_to_row.index:
                selected = locus_to_row.loc[override_locus]
                midpoint = (float(site["ChIP-exo Start"]) + float(site["ChIP-exo End"])) / 2.0
                computed_distance = abs(float(selected["tss"]) - midpoint)
                status = "mapped_curated_override"
                candidate_text = selected["gene_key"]
            else:
                status = "unresolved_override_locus_not_in_GFF"
        else:
            selected, status, candidate_text, computed_distance = choose_tu_gene(
                site, genes, gene_to_rows, gene_keys
            )

        row = {
            "Peak": site["Peak"],
            "Transcription_Unit": site["Transcription Unit"],
            "n_motifs": site["n_motifs"],
            "S_N_ratio": site["S/N ratio"],
            "paper_distance_to_tss": site.get("Distance to TSS", np.nan),
            "TU_gene_parse": candidate_text,
            "mapping_status": status,
            "matched_locus_tag": np.nan,
            "matched_gene": np.nan,
            "matched_start": np.nan,
            "matched_end": np.nan,
            "matched_strand": np.nan,
            "matched_tss": np.nan,
            "computed_distance_bp": computed_distance,
        }
        if selected is not None:
            row.update({
                "matched_locus_tag": selected["locus_tag"],
                "matched_gene": selected["gene"],
                "matched_start": selected["start"],
                "matched_end": selected["end"],
                "matched_strand": selected["strand"],
                "matched_tss": selected["tss"],
            })
        output_rows.append(row)

    out = pd.DataFrame(output_rows)
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    out.to_csv(OUT_PATH, sep="\t", index=False)

    mapped = out["mapping_status"].str.startswith("mapped_").sum()
    unresolved = len(out) - mapped
    print(f"Site-to-gene table written: {OUT_PATH}")
    print(f"Mapped: {mapped}/{len(out)}; unresolved: {unresolved}/{len(out)}")
    if unresolved:
        print("Review unresolved rows and add verified TU→locus_tag pairs to tu_primary_gene_map.tsv.")
        print(out.loc[~out["mapping_status"].str.startswith("mapped_"), [
            "Peak", "Transcription_Unit", "mapping_status", "TU_gene_parse"
        ]].to_string(index=False))


if __name__ == "__main__":
    main()
