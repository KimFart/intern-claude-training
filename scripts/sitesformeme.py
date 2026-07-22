import math
import os
import sys

import pandas as pd
from Bio import SeqIO

XLSX_PATH = "data/reference/ncomms5910_supplementary_data1.xlsx"
REFERENCE_DIR = "data/reference"
OUT_FASTA_PATH = "data/reference/fur_sites_for_meme.fasta"
FLANK = 50
FUR_BOX_CONSENSUS = "GATAATGATAATCATTATC"  # 19bp Fur box (Escolar et al. 1998, J Bacteriol)
FUR_BOX_HIT_THRESHOLD = 0.7


def resolve_reference_fasta(reference_dir=REFERENCE_DIR):
    """Return the path to the E. coli K-12 MG1655 reference FASTA in reference_dir,
    logging which RefSeq version was actually used.

    Seo et al. 2014's ChIP-exo coordinates were called against NC_000913.2.
    NC_000913.3 carries a handful of sequence corrections (indels) relative to
    .2 -- using .3 can shift coordinates near/downstream of any correction.
    Prefer .2 if present; warn loudly if falling back to .3, since that
    mismatch is exactly the kind of thing positive_control() below is meant to
    catch before it silently corrupts every downstream result.
    """
    candidates = ["NC_000913.2.fasta", "NC_000913.fasta", "NC_000913.3.fasta"]
    for name in candidates:
        path = os.path.join(reference_dir, name)
        if os.path.isfile(path):
            print(f"[genome] using reference FASTA: {path}")
            if "913.3" in name:
                print(
                    "[genome] WARNING: paper coordinates are on NC_000913.2; NC_000913.3 has "
                    "since-corrected indels that can shift coordinates. Check the "
                    "positive_control() result below before trusting any extraction."
                )
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


def _reverse_complement_str(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def best_fur_box_match(seq, consensus=FUR_BOX_CONSENSUS):
    """Slide the literal Fur-box consensus across seq (both strands) and
    return the single best-matching window as
    (identity_fraction, strand, window_start_offset, matched_subseq).

    This is a crude sanity check (simple positional identity, no IUPAC
    degeneracy or PWM) meant to catch gross problems -- wrong coordinates,
    wrong organism, sites that don't actually carry a Fur-box-like sequence --
    before spending time on MEME. It is not a substitute for MEME's own
    statistical motif model.
    """
    w = len(consensus)
    best = (0.0, "+", -1, "")
    for strand, candidate in (("+", seq), ("-", _reverse_complement_str(seq))):
        for i in range(len(candidate) - w + 1):
            window = candidate[i : i + w].upper()
            identity = sum(a == b for a, b in zip(window, consensus)) / w
            if identity > best[0]:
                best = (identity, strand, i, window)
    return best


def report_fur_box_hit_rate(df, genome_seq, strand_col, threshold=FUR_BOX_HIT_THRESHOLD):
    """Score every (not just top-S/N) extracted site against FUR_BOX_CONSENSUS
    and report how many clear `threshold` identity, plus the identity
    distribution -- so a real Fur-box shortfall is visible up front instead of
    only showing up later as a surprising MEME result.
    """
    scores = []
    for _, row in df.iterrows():
        strand = get_strand(row, strand_col)
        result = extract_site_sequence(genome_seq, row["ChIP-exo Start"], row["ChIP-exo End"], strand)
        if result is None:
            continue
        seq, _, _ = result
        identity, _, _, _ = best_fur_box_match(seq)
        scores.append(identity)

    hits = sum(s >= threshold for s in scores)
    scores.sort()
    print(
        f"\nFur-box consensus check ({FUR_BOX_CONSENSUS}, {len(FUR_BOX_CONSENSUS)}bp): "
        f"{hits}/{len(scores)} sites reach >={threshold:.0%} identity (best of both strands)"
    )
    print(
        f"Identity distribution: min={scores[0]:.0%}  median={scores[len(scores) // 2]:.0%}  "
        f"max={scores[-1]:.0%}"
    )


def positive_control(df, genome_seq, chrom_id, strand_col, n=5, max_mismatch=5,
                      box_col="Significance (p-value) of similarity to Fur box"):
    """Verify the extraction is actually correct before handing anything to MEME.

    Picks the n sites the paper itself scored as the strongest Fur-box matches
    (smallest p-value in box_col; falls back to a short list of textbook Fur
    targets if that column is unusable) and checks whether the extracted
    window actually contains something close to FUR_BOX_CONSENSUS
    (<= max_mismatch out of 19bp, either strand). If not even one of these
    known-strong sites passes, coordinates or genome version are wrong --
    extraction must be fixed before running MEME, so this exits hard rather
    than letting a silent bug flow into "surprising" motif results downstream.
    """
    if box_col in df.columns:
        targets = df.copy()
        targets[box_col] = pd.to_numeric(targets[box_col], errors="coerce")
        targets = targets.dropna(subset=[box_col]).sort_values(box_col).head(n)
    else:
        known = {"fepA-entD", "fhuACDB", "feoABC", "mntH", "fecIR"}
        targets = df[df["Transcription Unit"].isin(known)].head(n)

    print(f"\n--- Positive control: {len(targets)} known strong Fur-box sites ---")
    passed = 0
    for _, row in targets.iterrows():
        strand = get_strand(row, strand_col)
        result = extract_site_sequence(genome_seq, row["ChIP-exo Start"], row["ChIP-exo End"], strand)
        if result is None:
            print(f"  {row['Peak']}: out of chromosome bounds")
            continue
        seq, window_start, window_end = result
        identity, match_strand, _, match_seq = best_fur_box_match(seq)
        mismatches = round((1 - identity) * len(FUR_BOX_CONSENSUS))
        ok = mismatches <= max_mismatch
        passed += ok
        print(
            f"  {row['Peak']:5s} {str(row['Transcription Unit']):18s} "
            f"{chrom_id}:{window_start}-{window_end}({strand})  "
            f"best Fur-box = {mismatches}/{len(FUR_BOX_CONSENSUS)} mismatches  "
            f"strand={match_strand}  {match_seq}  -> {'PASS' if ok else 'FAIL'}"
        )
    print(f"--- {passed}/{len(targets)} passed (<= {max_mismatch} mismatches) ---")

    if passed == 0:
        sys.exit(
            "FATAL: none of the known strong Fur-box sites contain a recognizable Fur box.\n"
            "       Coordinates or genome version (NC_000913.2 vs .3) are likely mismatched.\n"
            "       Fix extraction before running MEME -- do not trust the output FASTA."
        )
    elif passed < len(targets):
        print(
            "NOTE: not all known-strong sites passed. Not a total failure, but worth "
            "checking the coordinates/window for the ones that failed."
        )


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

    report_fur_box_hit_rate(df, genome_seq, strand_col)
    positive_control(df, genome_seq, chrom_id, strand_col)


if __name__ == "__main__":
    main()
