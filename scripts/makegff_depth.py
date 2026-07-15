#!/usr/bin/env python3
# follow-up 1: column 1 now contains the element before '.' by splitting the reference name on '.' and taking the first part. This is to match the GFF annotation file which has the same format.
# follow-up 2: column 6 now put a '.' instead of hardcoding the score as 1.
# follow-up 3: column 6 now reports the number of reads sharing the same 5' position
#              and strand as this read (ChIP-exo border pileup depth), instead of '.'.
# follow-up 4: column 6 is now signed by strand (negative on '-'), so a single
#              track mirrors + depth above the axis and - depth below it,
#              matching the strand-paired plots in ChIP-exo papers. (Dropped
#              the local-maxima/find_peaks version from follow-up 4 - too much
#              complexity for this script; every read's row is scored now,
#              not just peaks.)
"""Convert a sorted BAM file to a GFF file of aligned reads (one row per read),
scoring each row by the number of reads sharing its 5' position and strand
(the ChIP-exo border pileup depth at that read's origin). Scores on the '-'
strand are negated so the two strands mirror around a zero baseline.
"""

import argparse
import os
import sys
from collections import Counter

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a sorted BAM file to a GFF file of aligned reads, "
        "scored by strand-signed 5'-position read depth."
    )
    parser.add_argument("bam_filepath", help="path to the input sorted .bam file")
    parser.add_argument("gff_filepath", help="path to write the output .gff file")
    return parser.parse_args()


def five_prime_pos(read):
    """0-based 5' position of a mapped read, accounting for strand."""
    return read.reference_end - 1 if read.is_reverse else read.reference_start


def count_five_prime_depth(bam):
    """Return {(seqname, strand, five_prime_pos): read_count} over all mapped reads."""
    counts = Counter()
    for read in bam:
        if read.is_unmapped:
            continue
        seqname = read.reference_name.split(".")[0]
        strand = "-" if read.is_reverse else "+"
        counts[(seqname, strand, five_prime_pos(read))] += 1
    bam.reset()
    return counts


def main():
    args = parse_args()

    if not os.path.isfile(args.bam_filepath):
        sys.exit(f"Error: file not found: {args.bam_filepath}")

    total_reads = 0
    written_reads = 0
    skipped_unmapped = 0

    with pysam.AlignmentFile(args.bam_filepath, "rb") as bam:
        depth_counts = count_five_prime_depth(bam)

        with open(args.gff_filepath, "w") as gff:
            for read in bam:
                total_reads += 1
                if read.is_unmapped:
                    skipped_unmapped += 1
                    continue

                seqname = read.reference_name.split(".")[0]
                start = read.reference_start + 1
                end = read.reference_end
                strand = "-" if read.is_reverse else "+"
                name = read.query_name
                depth = depth_counts[(seqname, strand, five_prime_pos(read))]
                score = -depth if strand == "-" else depth

                gff.write(f"{seqname}\tmakegff\tread\t{start}\t{end}\t{score}\t{strand}\t.\tname={name}\n")
                written_reads += 1

    print(f"Total reads read:         {total_reads}")
    print(f"Reads skipped (unmapped): {skipped_unmapped}")
    print(f"Reads written:            {written_reads}")
    print(f"Output written to:        {args.gff_filepath}")


if __name__ == "__main__":
    main()
