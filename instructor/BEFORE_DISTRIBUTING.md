# Instructor Checklist — Before Distributing to Interns

Complete every item below before sharing the Codespace link with a new cohort.

---

## 1. SRR Accession — Module 3

**File to edit:** `data/reference/README.md`

Replace both instances of `[SRR_ACCESSION]` with the actual accession number for your cohort's dataset.

Requirements for the chosen accession:
- Single-end ChIP-exo reads from *E. coli* K-12 MG1655
- ~500,000 reads (subsample if necessary with `seqtk sample`)
- Alignment should complete in under 5 minutes on Codespace free tier (2 vCPU, 4 GB RAM)

Recommended: use a single-end ChIP-exo sample from the same Seo et al. 2014 study (GEO: GSE54901) for consistency with Module 4. Subsample with *(run this on the lab server, not inside the Codespace — `seqtk` is not installed in the intern environment)*:
```bash
seqtk sample SRR######.fastq 500000 > SRR######_500k.fastq
```

---

## 2. Lab Annotation GFF — All Modules

**File to place:** `data/reference/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff`

This file is the lab's *E. coli* K-12 MG1655 annotation GFF with TSS positions. It is not in the repo. Copy it from the lab server before distributing:

```bash
cp /path/to/lab/data/ec_annotation_20100903_DHK_cSRNA_with_ortho.gff \
   data/reference/
```

Used in Module 2 (pandas exercises), Module 4 Exercise 7 (TSS distance analysis), and the mini-project.

---

## 3. Mini-Project Dataset

No action needed. The mini-project uses the same Fur ChIP-exo files interns download themselves in Module 4.

---

## 4. Verify GSE54901 Download Works

Before distributing, manually open the GEO series page and confirm the peak file is available:

```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54901
```

- Confirm the Fur ChIP-exo peak file for the iron-replete condition is listed under **Supplementary file** and is downloadable
- Note the exact filename here: `[filename]` — record this in `data/reference/README.md` so interns can verify they got the right file
- If the file format has changed (NCBI occasionally restructures GEO supplementary files), update the Module 4 Exercise 4 instructions accordingly

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
5. Run Module 1 notebook — concept cells should render without errors. **Note:** Exercise 2 (`m1-ex2-code`) is intentionally broken and will raise a `KeyError` when run — this is expected and is the exercise.
6. Confirm plan mode works: in the Claude Code terminal, press **Shift+Tab** before sending a prompt

---

## 6. Claude Code Pro Plan

Interns need a Claude.ai Pro plan account to use plan mode (Shift+Tab).

- Confirm each intern has a Pro plan account before Week 3 (Module 3 uses plan mode for the first time)
- Interns authenticate by running `claude` in the Codespace terminal and following the login prompt

---

## 7. Final Checklist

- [ ] `[SRR_ACCESSION]` replaced in `data/reference/README.md`
- [ ] `ec_annotation_20100903_DHK_cSRNA_with_ortho.gff` copied to `data/reference/`
- [ ] GSE54901 supplementary file verified downloadable; filename recorded in `data/reference/README.md`
- [ ] Codespace test-launched and all tools verified (including `pysam`)
- [ ] Module 1 concept cells render correctly; Exercise 2 intentionally raises `KeyError`
- [ ] Plan mode tested in Codespace terminal
- [ ] All interns have Claude Pro plan accounts
- [ ] `instructor/rubric.md` reviewed and shared with interns (or kept internal — your call)
