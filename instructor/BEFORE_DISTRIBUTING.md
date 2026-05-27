# Instructor Checklist — Before Distributing to Interns

Complete every item below before sharing the Codespace link with a new cohort.

---

## 1. SRR Accession — Module 3

**File to edit:** `module3-ngs-pipeline/data/sample/README.md`

Replace both instances of `[SRR_ACCESSION]` with the actual accession number for your cohort's dataset.

Requirements for the chosen accession:
- Single-end ChIP-exo reads from *E. coli* K-12 MG1655
- ~500,000 reads (subsample if necessary with `seqtk sample`)
- Alignment should complete in under 5 minutes on Codespace free tier (2 vCPU, 4 GB RAM)

Recommended: pick from published FUR ChIP-exo datasets (e.g., GEO series GSE41d698 or similar). Subsample with:
```bash
seqtk sample SRR######.fastq 500000 > SRR######_500k.fastq
```

---

## 2. `makegff.py` — Module 3, Question 3-6

**Status:** This script is a lab custom tool and is not included in this repo.

Before distributing, either:
- Copy `makegff.py` from the lab server into `module3-ngs-pipeline/` (recommended)
- Or update Step 6 in `03_ngs_pipeline_with_claude.ipynb` to use `bedtools bamtobed` or another accessible alternative

The script is called in the pipeline as:
```bash
python makegff.py input.bam output.gff
```

---

## 3. Module 4 Paper Excerpt — FUR Binding Sites

**File:** `module4-multiomics-motif/data/sample/fur_paper_excerpt.txt`

This file contains a synthetic excerpt describing FUR binding site coordinates. Verify that the coordinates in the excerpt are biologically reasonable and consistent with the published FUR regulon (e.g., Seo et al. 2014, Nucleic Acids Res).

If you swap in a real paper excerpt, update the validation expectations accordingly.

---

## 4. Mini-Project Dataset — Modules 5–6

**Directory:** `module5-6-miniproject/data/`

This directory is empty in the repo. Add the dataset before distributing:
- Recommended: a subset of a real RNA-seq or ChIP-seq experiment from the lab
- Include a brief `data/README.md` explaining the organism, condition, and experiment type (without giving away the biological answer)
- Do not include raw FASTQ files larger than ~100 MB — use processed files (e.g., count matrices, peak BED files) where possible

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
   meme --version
   python -c "import biopython; print('ok')"
   efetch -version
   claude --version
   ```
5. Run Module 1 notebook top-to-bottom — it should produce no errors
6. Confirm plan mode works: in the Claude Code terminal, press **Shift+Tab** before sending a prompt

---

## 6. Claude Code Pro Plan

Interns need a Claude.ai Pro plan account to use plan mode (Shift+Tab).

- Confirm each intern has a Pro plan account before Week 3 (Module 3 uses plan mode for the first time)
- Interns authenticate by running `claude` in the Codespace terminal and following the login prompt

---

## 7. Final Checklist

- [ ] `[SRR_ACCESSION]` replaced in `module3-ngs-pipeline/data/sample/README.md`
- [ ] `makegff.py` copied to `module3-ngs-pipeline/` (or Step 6 updated)
- [ ] `module5-6-miniproject/data/` contains the cohort dataset
- [ ] Codespace test-launched and all tools verified
- [ ] Plan mode tested in Codespace terminal
- [ ] All interns have Claude Pro plan accounts
- [ ] `instructor/rubric.md` reviewed and shared with interns (or kept internal — your call)
