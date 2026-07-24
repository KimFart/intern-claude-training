#!/usr/bin/env python
# makegff.py — SBML lab BAM->GFF coverage converter (provided script).
# Counts per-base read coverage in a strand-specific manner and writes one GFF row
# per non-zero position, so pile-ups render as peaks in MetaScope.
#   --separate_strand : + and - strands as separate tracks
#   --flip            : swap strand assignment (use for RNA-seq; omit for ChIP-exo)
#   --log_scale       : write log2(count) instead of raw depth
# Output seqname is hardcoded NC_000913 (matches the lab annotation).
# Source: github.com/SKo-sierrayk/sbml_kimslab_tutorials (one-line py3 filename-decode fix).
# Usage:  python makegff.py [--separate_strand] [--flip] [--log_scale] <sorted.bam> <out.gff>

from os.path import split

from numpy import zeros

import math
import pysam


def count_coverage(samfile, chromosome_size=6000000, flip=False):
    """counts coverage per base in a strand-specific manner

    For paired-end reads, the insert betwen the mapped reads is
    also counted.

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    chromsome_size: This value should be larger than the largest chromosome"""
    all_counts = {}
    plus_strands = []
    minus_strands = []
    if "SQ" in samfile.header:
        chromosome_sizes = {}
        for entry in samfile.header["SQ"]:
            chromosome_sizes[entry["SN"]] = int(entry["LN"]) + 1
    else:
        for reference in samfile.references:
            chromosome_sizes[reference] = chromosome_size
    for reference in samfile.references:  # create an array for each reference
        plus_strands.append(zeros((chromosome_sizes[reference],)))
        minus_strands.append(zeros((chromosome_sizes[reference],)))
    # iterate through each mapped read
    for i, read in enumerate(samfile):
        if read.is_unmapped:
            continue
        if not read.is_proper_pair:
            if read.is_reverse:
                minus_strands[read.tid][read.pos:read.aend] += 1
            else:
                plus_strands[read.tid][read.pos:read.aend] += 1
        # for paired-end data, get entire insert from only read1
        elif read.is_read1:
            if read.is_reverse:
                minus_strands[read.tid][read.pnext:read.aend] += 1
            else:
                plus_strands[read.tid][read.pos:read.pos + read.isize] += 1
    # store the results per reference
    for i, reference in enumerate(samfile.references):
        all_counts[reference] = {}
        if flip:
            all_counts[reference]["-"] = plus_strands[i]
            all_counts[reference]["+"] = minus_strands[i]
        else:
            all_counts[reference]["+"] = plus_strands[i]
            all_counts[reference]["-"] = minus_strands[i]
    return all_counts

def count_coverage_indexed(samfile, chromosome_size=6000000, flip=False):
    """counts coverage per base in a strand-specific manner

    For paired-end reads, the insert betwen the mapped reads is
    also counted.

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    chromsome_size: This value should be larger than the largest chromosome"""
    if "SQ" in samfile.header:
        chromosome_sizes = {}
        for entry in samfile.header["SQ"]:
            chromosome_sizes[entry["SN"]] = int(entry["LN"]) + 1
    else:
        for reference in samfile.references:
            chromosome_sizes[reference] = chromosome_size
    all_counts = {}
    for reference in samfile.references:  # go through each chromosome
        plus_strand = zeros((chromosome_sizes[reference],))
        minus_strand = zeros((chromosome_sizes[reference],))
        # iterate through each mapped read
        for i, read in enumerate(samfile.fetch(reference=reference)):
            if not read.is_proper_pair:
                if read.is_reverse:
                    minus_strand[read.pos:read.aend] += 1
                else:
                    plus_strand[read.pos:read.aend] += 1
            # for paired-end data, get entire insert from only read1
            elif read.is_read1:
                if read.is_reverse:
                    minus_strand[read.pnext:read.aend] += 1
                else:
                    plus_strand[read.pos:read.pos + read.isize] += 1
            all_counts[reference] = {}
            if flip:
                all_counts[reference]["-"] = plus_strand
                all_counts[reference]["+"] = minus_strand
            else:
                all_counts[reference]["+"] = plus_strand
                all_counts[reference]["-"] = minus_strand
    return all_counts


def write_samfile_to_gff(samfile, output, chromosome_size=6000000,
                         separate_strand=False, flip=False, log_scale=False):
    """write samfile object to an output object in a gff format

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo

    chromsome_size: This value should be larger than the largest chromosome

    separate_strand: Whether the forward and reverse strands should be made
    into separate tracks (True) or the negative strand should be rendered
    as negative values (False)
    """
    all_counts = count_coverage(samfile, chromosome_size=chromosome_size,
                                flip=flip)
    _fn = samfile.filename
    if isinstance(_fn, bytes):
        _fn = _fn.decode()
    name = split(_fn)[1]
    for reference in all_counts:
        for strand in all_counts[reference]:
            counts = all_counts[reference][strand]			
            for i in counts.nonzero()[0]:			
                if log_scale:
                    count = math.log(float(counts[i]), 2)
                else:
                    count = counts[i]
                if separate_strand:
                    output.write("%s\t\t%s\t%d\t%d\t%.2f\t%s\t.\t.\n" %
                        ("NC_000913", "%s_(%s)" % (name, strand), i, i,count, strand))
                else:
                    output.write("%s\t\t%s\t%d\t%d\t%.2f\t%s\t.\t.\n" %
                        ("NC_000913", name, i, i, count,strand))
						

def convert_samfile_to_gff(sam_filename, out_filename, chromosome_size=6000000,
                           separate_strand=False, flip=False, log_scale=False):
    """read in the a samfile from a path, and write it out to a gff filepath

    flip: Whether or not the strands should be flipped.
    This should be true for RNA-seq, and false for ChIP-exo.

    chromsome_size: This value should be larger than the largest chromosome.

    separate_strand: Whether the forward and reverse strands should be made
    into separate tracks (True) or the negative strand should be rendered
    as negative values (False)
    """
    samfile = pysam.Samfile(sam_filename)
    with open(out_filename, "w") as outfile:
        write_samfile_to_gff(samfile, outfile, chromosome_size=chromosome_size,
                             separate_strand=separate_strand, flip=flip, 
							 log_scale=log_scale)
    samfile.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    from os.path import isfile

    parser = ArgumentParser("convert a samfile to a gff file")

    parser.add_argument("sam_filename")
    parser.add_argument("out_filename", nargs="?", default=None)
    parser.add_argument("--chromosome_size", type=int, default=6000000,
        help="""This value should be larger than the largest chromosome.
        Default is 6000000""")
    parser.add_argument("--flip", required=False, action="store_true",
        help="""Whether or not the strands should be flipped.
        This should be true for RNA-seq, and false for ChIP-exo.
        False by default.""")
    parser.add_argument("--separate_strand", required=False,
        action="store_true", help="""Whether the forward and reverse strands
            should be made into separate tracks. Otherwise the negative strand
            will be rendered as negative values. Default is False.""") 
    parser.add_argument("--log_scale", required=False,
        action="store_true", help="""Whether the count number should be
            in log-scale. Default is False.""")
    args = parser.parse_args()
    if args.out_filename is None:
        if args.sam_filename.endswith(".sam") or \
                args.sam_filename.endswith(".bam"):
            new_filename = args.sam_filename[:-3] + "gff"
            if isfile(new_filename):  # do not want to overwrite existing
                raise IOError("File %s already exists" % new_filename)
            else:
                args.out_filename = new_filename
    convert_samfile_to_gff(args.sam_filename, args.out_filename,
        chromosome_size=args.chromosome_size,
        separate_strand=args.separate_strand, flip=args.flip, log_scale=args.log_scale)
