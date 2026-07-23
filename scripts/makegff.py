#!/usr/bin/env python3

# makegff.py 스크립트는 sorted BAM을 입력으로 받아, samtools view -F 0xD04로 unmapped/secondary/duplicate FLAG가 켜진 리드를 제거한 뒤, 
# bedtools genomecov -5 -bg를 +/- strand에 각각 실행하여 각 리드의 5' 말단이 겹치는 depth(coverage)를 bedGraph 형식으로 계산한다.
# 이 bedGraph 결과를 0-based에서 GFF의 1-based 좌표로 변환하고, seqname은 '.'을 기준으로 split하여 첫 번째 부분만 ㄴ남겨 annotation GFF의 seqname과 일치시킨다.
# -strand의 depth 값은 부호를 반전시켜 MetaScope에서 하나의 트랙 안에서 +/- 신호를 시각적으로 구분할 수 있도록 하였다.
# --cap-percentile 옵션으로 지정한 백분위수 이상의 depth 값을 잘라내 outlier에 의한 시각적 왜곡을 줄일 수 있고, 
# --log-transform 옵션으로 log2(depth+1) 변환을 적용해 강한 신호를 정규화할 수 있으며, 변환 후 값은 소수점 4자리로 반올림된다.
# 최종적으로 seqname, source(makegff), feature(fiveprime), start, end, score(depth), strand, frame(.), attribute(depth=<value>) 형식의 GFF를 출력한다.

# follow-up 1: column 1 now contains the element before '.' by splitting the reference name on '.' and taking the first part.
#              This is to match the GFF annotation file which has the same format.
# follow-up 2: column 6 now contains the read-depth (coverage) at that position, computed via
#              bedtools genomecov -5, instead of a hardcoded score or '.'.
#              Minus-strand depth values are negated so that +/- strand signal can be
#              distinguished within a single track (MetaScope supports negative scores).
# follow-up 3: only the 5' end of each read is represented (ChIP-exo signal), not the full aligned span.
#              Depth counting is delegated to bedtools (subprocess) rather than reimplemented in pure Python.
# follow-up 4: depth values can optionally be normalized to reduce the visual impact of outliers,
#              using percentile-based capping (winsorizing) and/or log2 transform, applied on the
#              absolute depth value (sign is preserved to keep the strand encoding intact).
#              log2-transformed values are rounded to 4 decimal places (hardcoded, not a CLI flag)
#              since raw depth is an integer count -- full float precision (e.g. 3.169925001442312)
#              is meaningless noise, not extra information.
# follow-up 5: --separate_strand disables the minus-strand sign negation (used by Module 3/ChIP-exo
#              to pack +/- signal into one score column). RNA-seq coverage (Module 4) needs true,
#              unsigned depth on both strands, distinguished only by the GFF strand column -- negating
#              minus-strand RNA-seq depth would misrepresent it as a trough instead of real coverage.

"""Convert a sorted BAM file to a GFF file of 5'-end read-depth signal (ChIP-exo style)."""

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a sorted BAM file to a GFF file of aligned reads."
    )
    parser.add_argument("bam_filepath", help="path to the input sorted .bam file")
    parser.add_argument("gff_filepath", help="path to write the output .gff file")
    parser.add_argument(
        "--cap-percentile", type=float, default=None,
        help="cap depth values at this percentile (e.g. 99) to reduce outlier impact; default: no capping",
    )
    parser.add_argument(
        "--log-transform", action="store_true",
        help="apply log2(depth + 1) transform after capping",
    )
    parser.add_argument(
        "--separate_strand", action="store_true",
        help="keep minus-strand depth positive instead of negating it, so +/- signal is "
             "distinguished by the GFF strand column rather than by score sign (use for RNA-seq)",
    )
    return parser.parse_args()


def check_bedtools():
    if shutil.which("bedtools") is None:
        sys.exit("Error: 'bedtools' executable not found in PATH. Please install bedtools.")


def check_samtools():
    if shutil.which("samtools") is None:
        sys.exit("Error: 'samtools' executable not found in PATH. Please install samtools.")


def run_genomecov(bam_filepath, strand, out_path):
    """Run bedtools genomecov -5 -bg -strand <strand> and write result to out_path.

    Reads are piped through `samtools view -F 0xD04` first to drop unmapped,
    secondary, and duplicate-flagged alignments -- bedtools genomecov -ibam
    does not filter these on its own and would otherwise inflate depth.
    """
    view_cmd = ["samtools", "view", "-b", "-F", "0xD04", bam_filepath]
    genomecov_cmd = [
        "bedtools", "genomecov",
        "-5",              # count only the 5' end of each read
        "-bg",              # bedGraph output (chrom, start, end, depth)
        "-strand", strand,
        "-ibam", "stdin",
    ]
    with open(out_path, "w") as out_f:
        view_proc = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
        subprocess.run(genomecov_cmd, stdin=view_proc.stdout, stdout=out_f, check=True)
        view_proc.stdout.close()
        if view_proc.wait() != 0:
            sys.exit(f"Error: samtools view failed on {bam_filepath}")


def bedgraph_to_gff_rows(bedgraph_path, strand_symbol, negate=False):
    """Yield (seqname, start, end, depth, strand) tuples from a bedGraph file."""
    with open(bedgraph_path) as f:
        for line in f:
            chrom, start, end, depth = line.rstrip("\n").split("\t")
            seqname = chrom.split(".")[0]
            gff_start = int(start) + 1   # bedGraph is 0-based half-open -> GFF is 1-based inclusive
            gff_end = int(end)
            depth_val = -int(depth) if negate else int(depth)
            yield (seqname, gff_start, gff_end, depth_val, strand_symbol)


def compute_percentile_cap(abs_depths, percentile):
    """Return the depth value at the given percentile of the (sorted) abs_depths list."""
    if not abs_depths:
        return 0
    sorted_depths = sorted(abs_depths)
    idx = int(len(sorted_depths) * percentile / 100)
    idx = min(idx, len(sorted_depths) - 1)
    return sorted_depths[idx]


LOG_TRANSFORM_ROUND_DECIMALS = 4


def normalize_depth(depth_val, cap=None, log_transform=False):
    """Apply optional percentile cap and/or log2 transform, preserving the original sign."""
    sign = -1 if depth_val < 0 else 1
    value = abs(depth_val)

    if cap is not None:
        value = min(value, cap)
    if log_transform:
        value = math.log2(value + 1)
        value = round(value, LOG_TRANSFORM_ROUND_DECIMALS)

    return sign * value

def main():
    args = parse_args()

    if not os.path.isfile(args.bam_filepath):
        sys.exit(f"Error: file not found: {args.bam_filepath}")

    check_bedtools()
    check_samtools()

    written_rows = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        plus_bg = os.path.join(tmpdir, "plus.bedgraph")
        minus_bg = os.path.join(tmpdir, "minus.bedgraph")

        run_genomecov(args.bam_filepath, "+", plus_bg)
        run_genomecov(args.bam_filepath, "-", minus_bg)

        rows = (
            list(bedgraph_to_gff_rows(plus_bg, "+", negate=False))
            + list(bedgraph_to_gff_rows(minus_bg, "-", negate=not args.separate_strand))
        )
        rows.sort(key=lambda r: (r[0], r[1]))  # sort by seqname, start

        cap = None
        if args.cap_percentile is not None:
            abs_depths = [abs(r[3]) for r in rows]
            cap = compute_percentile_cap(abs_depths, args.cap_percentile)
            print(f"Depth cap at {args.cap_percentile}th percentile: {cap}")

        with open(args.gff_filepath, "w") as gff:
            for seqname, start, end, depth, strand in rows:
                norm_depth = normalize_depth(depth, cap=cap, log_transform=args.log_transform)
                gff.write(f"{seqname}\tmakegff\tfiveprime\t{start}\t{end}\t{norm_depth}\t{strand}\t.\tdepth={norm_depth}\n")
                written_rows += 1

    print(f"GFF rows written:   {written_rows}")
    print(f"Output written to:  {args.gff_filepath}")


if __name__ == "__main__":
    main()
