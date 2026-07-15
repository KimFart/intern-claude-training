#!/usr/bin/env python3
# follow-up 1: column 1 now contains the element before '.' by splitting the reference name on '.' and taking the first part. This is to match the GFF annotation file which has the same format.
# follow-up 2: column 6 now put a '.' instead of hardcoding the score as 1.
"""Convert a sorted BAM file to a GFF file of aligned reads (one row per read)."""

import argparse
import os
import sys

import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a sorted BAM file to a GFF file of aligned reads."
    )
    parser.add_argument("bam_filepath", help="path to the input sorted .bam file")
    parser.add_argument("gff_filepath", help="path to write the output .gff file")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.bam_filepath):
        sys.exit(f"Error: file not found: {args.bam_filepath}")

    total_reads = 0
    written_reads = 0
    skipped_unmapped = 0

    with pysam.AlignmentFile(args.bam_filepath, "rb") as bam, open(args.gff_filepath, "w") as gff:
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

            gff.write(f"{seqname}\tmakegff\tread\t{start}\t{end}\t.\t{strand}\t.\tname={name}\n")
            written_reads += 1

    print(f"Total reads read:         {total_reads}")
    print(f"Reads skipped (unmapped): {skipped_unmapped}")
    print(f"Reads written:            {written_reads}")
    print(f"Output written to:        {args.gff_filepath}")


if __name__ == "__main__":
    main()
