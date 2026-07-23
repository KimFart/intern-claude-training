# follow-up 1: added timing around the iterrows() loop, printed at the end
# follow-up 2: also report the closest feature's TSS position, not just its distance
import time

import pandas as pd

GFF_PATH = "/workspaces/intern-claude-training/data/reference/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff"
TARGET = 1_000_000

gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# Count lines in the GFF file before loading with pandas
with open(GFF_PATH) as f:
    total_lines = sum(1 for _ in f)
print(f"Total lines in GFF file: {total_lines}")

df = pd.read_csv(GFF_PATH, sep="\t", header=None, names=gff_cols)
assert len(df) == total_lines, "line count / row count mismatch"

closest_gene_name = None
closest_distance = None
closest_tss = None
next_milestone = 10

start_time = time.perf_counter()
for idx, row in df.iterrows():
    tss = row["start"] if row["strand"] == "+" else row["end"]
    distance = abs(tss - TARGET)

    if closest_distance is None or distance < closest_distance:
        closest_distance = distance
        closest_tss = tss
        closest_gene_name = row["attributes"].split("gene=")[1].split(";")[0]

    rows_done = idx + 1
    pct_done = 100 * rows_done / total_lines
    if pct_done >= next_milestone:
        print(f"Progress: {next_milestone}% ({rows_done}/{total_lines} rows)")
        next_milestone += 10

elapsed = time.perf_counter() - start_time

print(f"Closest gene to position {TARGET}: {closest_gene_name}, TSS={closest_tss}, distance={closest_distance}")
print(f"Time spent on iterrows() loop: {elapsed:.4f} seconds")
