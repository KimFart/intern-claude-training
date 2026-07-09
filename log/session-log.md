## Session — 2026-07-08

### Done
- Built gff_filter.py, a CLI script that filters a GFF file by feature type, position range, strand, and attributes, with matching -exclude-* NOT filters, an -outputdir option, and coordinate validation.

### Broke / Struggled
- Non-numeric start/end coordinate values were silently passing through the position filters with no warning, until it was found and fixed.

### Learned
- To run an argparse-based script inline in a Jupyter cell, sys.argv must be manually overridden before calling main(), since the Jupyter kernel's own argv otherwise interferes.

---
## Session — 2026-07-08

### Done
- Loaded the GFF file with pandas and fixed Exercise 2's filepath bug.

### Broke / Struggled
- The KeyError: 'strand' error still shows up.

### Learned
- Using Pandas can significantly improve parsing GFF files.

---
## Session — 2026-07-09

### Done
- Built a GFF parsing/filtering python script.
- Built a GFF reading script using pandas' read_csv function.
- Fixed a malfunctioning read_csv script.
- Used the iterrows() function to process each row's column information.

### Broke / Struggled
- The markdown cell for Exercise 3a was placed beneath the cell for Exercise 3b.

### Learned
- Learned what a GFF format is, how to use it, and how to interpret some of the information in it.
- Now understand how iterrows() works, when to use it, and how to utilize it using Claude Code.

---
