import math
import os
import sys

import pandas as pd
from Bio import SeqIO

XLSX_PATH = "data/reference/ncomms5910_supplementary_data1.xlsx"
REFERENCE_DIR = "data/reference"
OUT_FASTA_PATH = "data/reference/fur_sites_for_meme.fasta"
FLANK = 50


def resolve_reference_fasta(reference_dir=REFERENCE_DIR):
    """Return the path to the E. coli K-12 MG1655 reference FASTA in reference_dir.

    Tries both filenames seen across this repo/notebook (with and without the
    RefSeq version suffix) so the script doesn't break depending on how it was
    downloaded.
    """
    candidates = ["NC_000913.3.fasta", "NC_000913.fasta"]
    for name in candidates:
        path = os.path.join(reference_dir, name)
        if os.path.isfile(path):
            return path
    tried = [os.path.join(reference_dir, n) for n in candidates]
    sys.exit(f"Error: reference FASTA not found. Tried: {tried}")


def load_binding_sites(xlsx_path=XLSX_PATH):
    """Load Seo et al. 2014 Supplementary Data 1 and return the iron-replete subset.

    The sheet has a title row above the real header, hence header=1. One row is
    an all-NaN placeholder (dropped via dropna(subset=["Peak"])). Iron-replete
    sites are rows where 'Binding Conditiona' contains 'R' (labels 'R' or 'R/S').
    """
    if not os.path.isfile(xlsx_path):
        sys.exit(f"Error: binding-site table not found: {xlsx_path}")
    df = pd.read_excel(xlsx_path, header=1)
    df = df.dropna(subset=["Peak"])
    print(f"Binding sites in table (after dropping placeholder row): {len(df)}")
    df = df[df["Binding Conditiona"].astype(str).str.contains("R")]
    print(f"Iron-replete sites (Binding Conditiona contains 'R'): {len(df)}")
    return df.reset_index(drop=True)


def find_strand_column(df):
    """Return the strand column name if the table has one, else None.

    Only literal 'Strand'/'strand' columns count. 'Binding Conditiona'
    ('R'/'S'/'R/S') encodes iron condition, not strand -- despite superficially
    resembling +/- labels, it must not be treated as one.
    """
    for col in ("Strand", "strand"):
        if col in df.columns:
            return col
    return None


def get_strand(row, strand_col):
    """Return '+' or '-' for a binding-site row.

    Falls back to '+' when the table has no strand column at all, or when a
    given row's strand value is missing.
    """
    if strand_col is None or pd.isna(row[strand_col]):
        return "+"
    return "-" if str(row[strand_col]).strip() == "-" else "+"


def site_window(start_1based, end_1based, flank=FLANK):
    """Compute the 0-based Python slice bounds for a +/-flank window around
    the midpoint of a 1-based inclusive [start, end] interval.

    Returns (slice_start, slice_end, window_start_1based, window_end_1based).
    """
    start_1based, end_1based = int(start_1based), int(end_1based)
    center = math.floor((start_1based + end_1based) / 2 + 0.5)  # explicit round-half-up
    window_start_1based = center - flank
    window_end_1based = center + flank
    # 1-based inclusive position p -> 0-based index (p - 1).
    # Slice start needs the -1; slice end doesn't, because Python's exclusive
    # slice end already lines up with the 1-based inclusive end after the
    # -1 (to 0-based) / +1 (exclusive-end) cancellation.
    slice_start = window_start_1based - 1
    slice_end = window_end_1based
    return slice_start, slice_end, window_start_1based, window_end_1based


def extract_site_sequence(genome_seq, start_1based, end_1based, strand, flank=FLANK):
    """Slice the +/-flank window around a binding site out of genome_seq.

    Returns (sequence_str, window_start_1based, window_end_1based), reverse-
    complemented if strand == '-'. Returns None (rather than silently wrapping
    on a negative index) if the window runs off either end of the chromosome.
    """
    slice_start, slice_end, window_start, window_end = site_window(start_1based, end_1based, flank)
    if slice_start < 0 or slice_end > len(genome_seq):
        return None
    window_seq = genome_seq[slice_start:slice_end]
    if strand == "-":
        window_seq = window_seq.reverse_complement()
    return str(window_seq), window_start, window_end


def build_fasta_header(peak_id, chrom_id, window_start, window_end, strand):
    """Build a MEME-safe, traceable FASTA header for one extracted site.

    Uses the actual extracted window coordinates (not the raw ChIP-exo
    Start/End) so the record is self-describing even without the source table.
    No spaces, since FASTA/MEME parsers split headers on the first whitespace.
    """
    return f"{peak_id}_{chrom_id}:{window_start}-{window_end}({strand})"


def write_fasta(records, out_path):
    """Write (header, sequence) tuples to out_path as flat, unwrapped FASTA."""
    with open(out_path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")
    print(f"FASTA records written: {len(records)}")
    print(f"Output written to: {out_path}")


def sanity_check(df, genome_seq, chrom_id, strand_col, n=3):
    """Print the window sequence + %AT content for the n highest-S/N-ratio
    sites, so a human can eyeball AT-richness (Fur boxes are AT-rich) before
    trusting the full batch. Does not try to match a literal Fur-box consensus
    -- confirming the real motif is MEME's job, next.
    """
    top = df.sort_values("S/N ratio", ascending=False).head(n)
    print(f"\n--- Sanity check: top {n} sites by S/N ratio ---")
    for _, row in top.iterrows():
        strand = get_strand(row, strand_col)
        result = extract_site_sequence(genome_seq, row["ChIP-exo Start"], row["ChIP-exo End"], strand)
        if result is None:
            print(f"{row['Peak']}: window out of chromosome bounds, skipped")
            continue
        seq, window_start, window_end = result
        at_frac = sum(b in "ATat" for b in seq) / len(seq)
        print(
            f"{row['Peak']}  {chrom_id}:{window_start}-{window_end}({strand})  "
            f"S/N={row['S/N ratio']:.2f}  AT content={at_frac:.0%}"
        )
        print(f"  {seq}")


def main():
    df = load_binding_sites()
    strand_col = find_strand_column(df)
    if strand_col is None:
        print(
            "No strand column found in the binding-site table -- defaulting all "
            "sites to '+'. This does not bias downstream results: MEME is run "
            "with -revcomp in the next section, which searches both strands "
            "regardless of how the input FASTA strand was recorded."
        )

    fasta_path = resolve_reference_fasta()
    genome_record = next(SeqIO.parse(fasta_path, "fasta"))
    genome_seq = genome_record.seq
    chrom_id = genome_record.id

    records = []
    skipped = 0
    for _, row in df.iterrows():
        strand = get_strand(row, strand_col)
        result = extract_site_sequence(genome_seq, row["ChIP-exo Start"], row["ChIP-exo End"], strand)
        if result is None:
            skipped += 1
            print(f"Skipped {row['Peak']}: window out of chromosome bounds")
            continue
        seq, window_start, window_end = result
        header = build_fasta_header(row["Peak"], chrom_id, window_start, window_end, strand)
        records.append((header, seq))

    write_fasta(records, OUT_FASTA_PATH)
    print(f"Sites skipped (window out of chromosome bounds): {skipped}")

    sanity_check(df, genome_seq, chrom_id, strand_col)


if __name__ == "__main__":
    main()
