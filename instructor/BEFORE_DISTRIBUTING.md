# Instructor Checklist — Before Distributing to Interns

Complete every item below before sharing the Codespace link with a new cohort.

---

## 1. SRR Accession — Module 3

No action needed. Interns find the SRR accession themselves in Module 3 Exercise 4 by navigating to GEO series GSE54901 and identifying a single-end ChIP-exo sample for the iron-replete condition. They then run `fastq-dump` themselves as part of the exercise.

---

## 2. Lab Annotation GFF — All Modules

**File to place:** `data/reference/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff`

This file is the lab's *E. coli* K-12 MG1655 gene annotation GFF (all rows are `gene` features; a TSS is derived from each gene's start/strand, not stored as a separate feature). It is not in the repo. Copy it from the lab server before distributing:

```bash
cp /path/to/lab/data/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff \
   data/reference/
```

Used in Module 2 (pandas exercises), Module 4 Exercise 7 (TSS distance analysis), and the mini-project.

---

## 3. Mini-Project Dataset

No action needed. The mini-project uses the same Fur ChIP-exo files interns download themselves in Module 4.

---

## 4. Verify the FUR data sources Work

Two distinct data sources are used in Module 4/5 — **GEO holds the raw ChIP-exo signal (coverage tracks), NOT the discrete binding sites.** The binding-site coordinates come from the paper's Supplementary Data Excel. Verify both:

**(a) The binding-site table (Module 4 Ex4 input — the critical one):**
```
https://www.nature.com/articles/ncomms5910
```
- Confirm the **Supplementary Data** Excel that lists Fur ChIP-exo binding sites (columns `Transcription Unit`, `Peak`, `ChIP-exo Start`, `ChIP-exo End`, `Binding Condition`, `S/N ratio`, `Significance`, `Distance to TSS`; 144 sites, header on the 2nd row) is downloadable. Interns download it and read it with `pd.read_excel(..., header=1)`.
- If the publisher restructures the supplementary files, update Module 4 Exercise 4 accordingly. The devcontainer includes `openpyxl` so pandas can read `.xlsx`; confirm `import openpyxl` works in the `sbml` env.

**(b) The GEO series (Module 3 alignment context / raw signal):**
```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54901
```
- Confirm the series is reachable. Note the supplementary GFFs here are **per-base coverage tracks** (millions of rows), not peak lists — this is expected and is what Module 3 aligns/visualizes.

---

## 5. GitHub Codespace — Test Before Distributing

1. Open the repo on GitHub
2. Click **Code → Codespaces → New codespace on main**
3. Wait for the `.devcontainer/setup.sh` to finish
4. In the terminal, verify:
   ```bash
   conda activate sbml
   bowtie2 --version
   samtools --version
   meme -version
   python -c "import Bio; print('Biopython', Bio.__version__)"
   python -c "import pysam; print('pysam', pysam.__version__)"
   efetch -version
   claude --version
   ```
5. Run Module 1 notebook — concept cells should render without errors. **Note:** Exercise 2 (`m1-ex2-code`) is intentionally broken — the countdown prints `0 1 2 3 4` instead of `5 4 3 2 1` (wrong `range()` arguments). No exception is raised; the wrong output is the bug interns must find.
6. Confirm plan mode works: in the Claude Code terminal, press **Shift+Tab** (cycle it until the status line shows plan mode is on) before sending a prompt. Verify the exact key sequence and indicator against the current Claude Code docs when testing.

---

## 6. Claude Code Pro Plan

Interns need a paid Claude plan account to run Claude Code at all (plan mode itself is a general interaction mode, not a plan-gated feature).

- Confirm each intern has a paid Claude plan account before Week 3 (Module 3 is the first module that uses Claude Code heavily, in plan mode)
- Interns authenticate by running `claude` in the Codespace terminal and following the login prompt

---

## 7. Final Checklist

- [ ] `ec_annotation_20100903_DHK_cSRNA_with_ortho.gff` copied to `data/reference/`
- [ ] FUR binding-site Excel verified downloadable from the Nature Communications article (Module 4 Ex4 input)
- [ ] GEO series GSE54901 reachable (raw ChIP-exo coverage signal — Module 3)
- [ ] GEO series GSE54900 reachable with a downloadable paired-end RNA-seq run (Module 4 §5 — e.g. SRR1168133 gives `_1`/`_2` FASTQ)
- [ ] `openpyxl` present in the `sbml` env (`pd.read_excel` on `.xlsx`)
- [ ] Codespace test-launched and all tools verified (including `pysam`)
- [ ] Module 1 concept cells render correctly; Exercise 2 intentionally prints wrong countdown output (not an exception)
- [ ] Plan mode tested in Codespace terminal
- [ ] All interns have Claude Pro plan accounts
- [ ] `instructor/rubric.md` reviewed and shared with interns (or kept internal — your call)
