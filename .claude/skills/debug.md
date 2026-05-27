---
name: debug
description: Structured error recovery. Use whenever a script or command fails. Enforces reproduce → localize → fix → document, preventing random guessing.
---

# Structured Debugging

## Process

Work through these steps in order. Do not skip to a fix before completing the diagnosis.

### Step 1 — Reproduce

Ask the user to provide:
- The exact command or code that failed
- The exact error message (full traceback if Python, full stderr if bash)
- What they expected to happen vs what actually happened

If anything is missing, ask for it before proceeding.

### Step 2 — Localize

Identify the root cause:
- What line or command is the actual failure point?
- Is it a file path issue, a format mismatch, a missing dependency, a wrong parameter, or a logic error?
- State your hypothesis clearly: "I think this is failing because..."

### Step 3 — Fix

Propose the minimal fix:
- Change only what is necessary
- Do not refactor or clean up surrounding code
- Show the corrected command or code snippet

### Step 4 — Verify

Tell the user how to confirm the fix worked:
- What output or file should they see if it succeeded?
- What command can they run to verify?

### Step 5 — Document

Remind the user: "Add this to your `/log` entry under 'Broke / Struggled' — note what caused it and what fixed it."

## Notes

- Common failure patterns in this lab's stack:
  - bowtie2: index not built before alignment (`bowtie2-build` must run first)
  - MEME: input FASTA must have ≥2 sequences; sequences shorter than `-maxw` will error — reduce `-maxw` or extend the extraction window in your Biopython script
  - Biopython SeqIO: wrong format string (e.g. `"fasta"` vs `"fastq"`)
  - fastq-dump: `--split-files` needed for paired-end; omitting it merges both reads into one file
