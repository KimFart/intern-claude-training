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
## Session — 2026-07-13

### Done
- Built scripts/run_pipeline.sh and scripts/makegff.py end-to-end (bowtie2-build → align → samtools sort/index → GFF), verified against synthetic and real data.
- Filled in notebook answers for the FASTQ/FASTA/SAM/BAM/GFF format question and Exercises 1-2.

### Broke / Struggled
- A real pipeline run against NC_000913.3 + ERR1254532 was killed mid-alignment with exit code 143 (SIGTERM) for an unclear reason.
- An efetch download printed a scary SSL/curl error but had actually succeeded via automatic retry.

### Learned
- Understood 0-based vs 1-based coordinate conversion (pysam reference_start/reference_end), why BAM sort must precede index, SAM FLAG bits, and FASTQ Phred quality scores.

---
## Session — 2026-07-14

### Done
- Completed Module 3 Exercises 1-7 end-to-end.
- Wrote and verified scripts/makegff.py (BAM to GFF, matches spec exactly).
- Built scripts/download_chipexo_reads.sh and scripts/download_reference_genome.sh.
- Added thread-count flags to scripts/run_pipeline.sh and the notebook's Exercise 7 cell.
- Ran the real pipeline (bowtie2 -> samtools -> makegff.py), producing a verified SRR1168122_chipexo.gff (685,127 reads, 99.22% alignment rate) with matching output in both the terminal and the notebook kernel.

### Broke / Struggled
- A notebook %%bash cell resolved bare `python` to the base conda env (no pysam) instead of the sbml env.

### Learned
- Many bash command files can be built almost automatically with Claude Code.
- %%writefile must be the first line of the cell.
- How pysam perceives a BAM record's FLAG.

---
## Session — 2026-07-15

### Done
- Nearly completed Module 3 and finished Exercise 10.
- Fixed makegff.py so its output GFF matches the annotation GFF's format and displays strand-sensitively in MetaScope.

### Broke / Struggled
- The Codespace build kept failing, requiring repeated Dockerfile edits and rebuilds; initially suspected a MEME Suite compatibility issue, but it turned out sratoolkit was failing to install (the rest of Module 3 wasn't blocked, so proceeded anyway).
- makegff.py wrote the versioned reference name NC_000913.3 (inherited from the reference FASTA header) while the annotation track used the unversioned NC_000913, silently splitting the two tracks onto separate MetaScope chromosome tabs instead of erroring.
- The GFF score column was hardcoded to 1; changed to . to match the annotation GFF and correctly separate strands in the MetaScope display.

### Learned
- Understood the principles of ChIP-exo sequencing in more depth.

---
## Session — 2026-07-15

### Done
- Reviewed makegff.py, then built scripts/makegff_depth.py (strand-signed 5'-end depth scoring) and scripts/run_pipeline_depth.sh, verified against the real SRR1168122 BAM.

### Broke / Struggled
- The scipy.signal.find_peaks-based local-maxima version worked correctly but was reverted for being too complicated for the script's purpose.

### Learned
- 5'-end (tag) counting vs full pileup depth for ChIP-exo border detection, and strand-signed mirrored score tracks as a common visualization convention.

---
## Session — 2026-07-16

### Done
- Rewrote makegff.py to delegate depth counting to bedtools genomecov instead of pure-Python pysam counting, and consolidated it into the canonical scripts/makegff.py + scripts/run_pipeline.sh pipeline.
- Fixed a real bug where bedtools genomecov was silently counting secondary and PCR-duplicate-flagged reads as real signal, inflating ChIP-exo depth; fixed by piping through samtools view -F 0xD04 first.
- Matched the GFF output format to the fixed notebook's Exercise 6 spec (fiveprime feature type, depth=<count> attribute) and renamed the output convention to SRR1168122_chipexo.gff.
- Fixed several notebook bugs in 03_ngs_pipeline_with_claude_fixed.ipynb (cells 22, 24, 29) caused by wrong assumptions about the kernel's working directory, and produced the Exercise 9/10 MetaScope deliverable PNG.

### Broke / Struggled
- Not related to the notebook itself: found several false-positive peaks that resemble an actual true-positive Fur peak; needs more digging into ChIP-exo literature to identify the data correctly.

### Learned
- How actual ChIP-exo pipeline elements quantify and normalize ChIP-exo coverage.

---
## Session — 2026-07-20

### Done
- Ran SRR1168133 (paired-end RNA-seq) through the pipeline: downloaded, aligned, sorted, indexed, and converted to GFF.
- Made slight modifications in makegff.py.

### Broke / Struggled
- Nothing broke — everything was okay today.

### Learned
- Why Fur is important to bacteria, and how it works in three different modes.

---
