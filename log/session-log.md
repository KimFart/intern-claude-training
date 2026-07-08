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
