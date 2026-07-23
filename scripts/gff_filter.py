#!/usr/bin/env python3

# Follow-up 1: added NOT logic for each filter type, so you can exclude features by type, position range, strand, and/or attributes.
# Follow-up 2: added a check for non-numeric start/end coordinates when using position filters, and a warning if they are present.i

"""Filter a 9-column GFF file by feature type, position range, strand, and/or attributes."""

import argparse
import os
import re
import sys

NUM_COLUMNS = 9
SEQID, SOURCE, TYPE, START, END, SCORE, STRAND, PHASE, ATTRIBUTES = range(NUM_COLUMNS)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter a GFF file by feature type, position range, strand, and/or attributes."
    )
    parser.add_argument("gff_filepath", help="path to the input .gff file")
    parser.add_argument("-featuretype", help='feature type to keep, e.g. "gene" (column 3)')
    parser.add_argument(
        "-positionrange",
        help='keep features with start > MIN and end < MAX, e.g. "500000-1000000" (columns 4/5)',
    )
    parser.add_argument("-strand", choices=["+", "-"], help="strand to keep (column 7)")
    parser.add_argument(
        "-outputdir",
        help="directory to write the filtered GFF to (default: same directory as the input file)",
    )
    parser.add_argument(
        "-attribute",
        action="append",
        dest="attributes",
        metavar="KEY=VALUE|KEY~SUBSTRING|KEY",
        help=(
            "attribute filter on column 9, repeatable (AND-combined). "
            'KEY=VALUE for exact match, KEY~SUBSTRING for substring match, '
            "KEY alone to require the key to be present. "
            'Examples: -attribute "gene=thrA" -attribute "product~kinase" -attribute "KpOrtho"'
        ),
    )
    parser.add_argument("-exclude-featuretype", help="feature type to drop, e.g. \"tRNA\" (column 3)")
    parser.add_argument(
        "-exclude-positionrange",
        help='drop features with start > MIN and end < MAX, e.g. "500000-1000000" (columns 4/5)',
    )
    parser.add_argument("-exclude-strand", choices=["+", "-"], help="strand to drop (column 7)")
    parser.add_argument(
        "-exclude-attribute",
        action="append",
        dest="exclude_attributes",
        metavar="KEY=VALUE|KEY~SUBSTRING|KEY",
        help="attribute filter on column 9 to drop, repeatable (AND-combined). Same syntax as -attribute.",
    )
    return parser.parse_args()


def parse_positionrange(value):
    parts = value.split("-")
    if len(parts) != 2:
        sys.exit(f"Error: -positionrange must look like 'MIN-MAX', got '{value}'")
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        sys.exit(f"Error: -positionrange must contain integers, got '{value}'")


def parse_attribute_filter(raw):
    if "~" in raw:
        key, _, value = raw.partition("~")
        return key.strip(), "~", value.strip()
    if "=" in raw:
        key, _, value = raw.partition("=")
        return key.strip(), "=", value.strip()
    return raw.strip(), "exists", None


def parse_attribute_column(text):
    attrs = {}
    for entry in text.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        key, sep, value = entry.partition("=")
        attrs[key.strip().lower()] = value.strip() if sep else ""
    return attrs


def attribute_matches(attrs_lower, key, op, value):
    key_lower = key.lower()
    if key_lower not in attrs_lower:
        return False
    if op == "exists":
        return True
    actual = attrs_lower[key_lower]
    if op == "=":
        return actual.lower() == value.lower()
    if op == "~":
        return value.lower() in actual.lower()
    return False


def sanitize(text):
    return re.sub(r"[^A-Za-z0-9+\-]", "_", text)


def attribute_token(key, op, value):
    return f"attr-{key}{op if op != 'exists' else ''}{value or ''}"


def build_output_path(
    input_path,
    featuretype,
    positionrange,
    strand,
    attribute_filters,
    exclude_featuretype,
    exclude_positionrange,
    exclude_strand,
    exclude_attribute_filters,
    outputdir=None,
):
    directory, filename = os.path.split(input_path)
    if outputdir:
        directory = outputdir
    stem, ext = os.path.splitext(filename)

    parts = []
    if featuretype:
        parts.append(sanitize(featuretype))
    if positionrange:
        start, end = positionrange
        parts.append(f"{start}-{end}")
    if strand:
        parts.append(strand)
    for key, op, value in attribute_filters:
        parts.append(sanitize(attribute_token(key, op, value)))
    if exclude_featuretype:
        parts.append(sanitize(f"not-{exclude_featuretype}"))
    if exclude_positionrange:
        start, end = exclude_positionrange
        parts.append(f"not-{start}-{end}")
    if exclude_strand:
        parts.append(f"not-{exclude_strand}")
    for key, op, value in exclude_attribute_filters:
        parts.append(sanitize(f"not-{attribute_token(key, op, value)}"))

    suffix = "_".join(parts)
    return os.path.join(directory, f"{stem}.{suffix}{ext}")


def parse_coords(fields):
    try:
        return int(fields[START]), int(fields[END])
    except ValueError:
        return None


def line_matches(
    fields,
    coords,
    featuretype,
    positionrange,
    strand,
    attribute_filters,
    exclude_featuretype,
    exclude_positionrange,
    exclude_strand,
    exclude_attribute_filters,
):
    if featuretype and fields[TYPE].lower() != featuretype.lower():
        return False
    if positionrange:
        if coords is None:
            return False
        start, end = coords
        min_pos, max_pos = positionrange
        if not (start > min_pos and end < max_pos):
            return False
    if strand and fields[STRAND] != strand:
        return False
    if attribute_filters:
        attrs_lower = parse_attribute_column(fields[ATTRIBUTES])
        for key, op, value in attribute_filters:
            if not attribute_matches(attrs_lower, key, op, value):
                return False

    if exclude_featuretype and fields[TYPE].lower() == exclude_featuretype.lower():
        return False
    if exclude_positionrange and coords is not None:
        start, end = coords
        min_pos, max_pos = exclude_positionrange
        if start > min_pos and end < max_pos:
            return False
    if exclude_strand and fields[STRAND] == exclude_strand:
        return False
    if exclude_attribute_filters:
        attrs_lower = parse_attribute_column(fields[ATTRIBUTES])
        for key, op, value in exclude_attribute_filters:
            if attribute_matches(attrs_lower, key, op, value):
                return False
    return True


def main():
    args = parse_args()

    if not os.path.isfile(args.gff_filepath):
        sys.exit(f"Error: file not found: {args.gff_filepath}")

    if args.outputdir and not os.path.isdir(args.outputdir):
        sys.exit(f"Error: output directory not found: {args.outputdir}")

    positionrange = parse_positionrange(args.positionrange) if args.positionrange else None
    attribute_filters = (
        [parse_attribute_filter(raw) for raw in args.attributes] if args.attributes else []
    )
    exclude_positionrange = (
        parse_positionrange(args.exclude_positionrange) if args.exclude_positionrange else None
    )
    exclude_attribute_filters = (
        [parse_attribute_filter(raw) for raw in args.exclude_attributes]
        if args.exclude_attributes
        else []
    )

    any_filter = bool(
        args.featuretype
        or positionrange
        or args.strand
        or attribute_filters
        or args.exclude_featuretype
        or exclude_positionrange
        or args.exclude_strand
        or exclude_attribute_filters
    )
    if not any_filter:
        print("No filter option added.")
        return

    needs_coords = bool(positionrange or exclude_positionrange)
    total_lines = 0
    skipped_lines = 0
    invalid_coord_lines = 0
    matched_lines = []

    with open(args.gff_filepath, newline="") as infile:
        for lineno, raw_line in enumerate(infile, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            total_lines += 1
            fields = line.split("\t")
            if len(fields) != NUM_COLUMNS:
                print(f"Warning: line {lineno} does not have {NUM_COLUMNS} columns, skipping")
                skipped_lines += 1
                continue
            coords = parse_coords(fields) if needs_coords else None
            if needs_coords and coords is None:
                print(f"Warning: line {lineno} has non-numeric start/end, position filters ignored for this line")
                invalid_coord_lines += 1
            if line_matches(
                fields,
                coords,
                args.featuretype,
                positionrange,
                args.strand,
                attribute_filters,
                args.exclude_featuretype,
                exclude_positionrange,
                args.exclude_strand,
                exclude_attribute_filters,
            ):
                matched_lines.append(line)

    output_path = build_output_path(
        args.gff_filepath,
        args.featuretype,
        positionrange,
        args.strand,
        attribute_filters,
        args.exclude_featuretype,
        exclude_positionrange,
        args.exclude_strand,
        exclude_attribute_filters,
        args.outputdir,
    )

    if os.path.exists(output_path):
        answer = input(f"'{output_path}' already exists. Overwrite? [y/N] ").strip().lower()
        if answer != "y":
            print("Aborted, no file written.")
            return

    with open(output_path, "w") as outfile:
        for line in matched_lines:
            outfile.write(line + "\n")

    print(f"Total lines read:      {total_lines}")
    print(f"Lines skipped:         {skipped_lines}")
    if needs_coords:
        print(f"Lines w/ bad coords:   {invalid_coord_lines}")
    print(f"Lines matched:         {len(matched_lines)}")
    print(f"Output written to:     {output_path}")


if __name__ == "__main__":
    main()
