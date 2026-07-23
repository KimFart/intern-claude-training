# SBML Intern Claude Code Training

> **Step 0 — Fork this repo first.** Before you do anything else, **fork this repository to your own GitHub account** and work in your fork — that is how you save and hand back your work (see [Submitting Your Work](#submitting-your-work)). Launch the Codespace **from your fork**, not from this template, or you won't be able to push your commits.

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/sbml-lab/intern-claude-training)

> The Codespace will install all required tools automatically (bowtie2, samtools, MEME, Biopython, pandas/openpyxl, Claude Code CLI). First launch takes ~5 minutes.
>
> You will need to authenticate Claude Code with a **Pro plan** account after launch: run `claude` in the terminal.

## Overview

A 6-week training curriculum for SBML Lab (KAIST GSEB) interns. Covers Python, GFF/pandas data processing, NGS alignment, motif analysis, and an independent mini-project — all using Claude Code as the primary tool.

Modules 3 and 4 close the loop by visualizing their GFF output in **MetaScope**, the lab's own genome browser. MetaScope is a desktop application (install it from the lab homepage, [sbml-lab.ai](https://sbml-lab.ai)) — the pipelines produce the GFF in the Codespace, and interns download it and open it in MetaScope locally to inspect ChIP-exo binding and RNA-seq signal against the gene annotation.

## Curriculum

| Module | Topic |
|--------|-------|
| 1 | What is Claude Code |
| 2 | GFF parsing + pandas with Claude Code |
| 3 | NGS alignment pipeline with Claude Code + MetaScope visualization of ChIP-exo |
| 4 | The FUR regulon: paper → paired-end RNA-seq alignment → MEME motif with Biopython → integrative MetaScope view |
| 5 (Weeks 5–6) | Independent mini-project |

## Standing Rules

All interns follow these rules throughout the training:

1. `/log` at the end of every session
2. `/debug` before asking for help when something breaks
3. `/explain [concept]` before Googling
4. Plan mode (Shift+Tab twice) before writing any pipeline

## Submitting Your Work

You hand back your work through **git**, not by sharing your Codespace. At the start, **fork this repository to your own GitHub account** and do all your work in your fork. When you're done (and at the end of every session), commit and push, then give your instructor **the URL of your fork** (`github.com/<your-username>/intern-claude-training`).

**What to commit — and where it must go:**

| Deliverable | Path (must be here to reach your instructor) |
|---|---|
| Completed notebooks | `notebooks/0*.ipynb` |
| Session logs (from `/log`) | `log/session-log.md` |
| Scripts you wrote (e.g. `makegff.py`) | `notebooks/` |
| **MetaScope exported figures** | **`notebooks/`** — e.g. `notebooks/module3_chipexo_metascope.png` |

> **⚠️ Two things that will silently NOT reach your instructor:**
> 1. **Anything in `data/reference/`** — raw reads, genomes, BAMs, and generated GFFs there are gitignored (too large to commit). That's fine; your instructor regenerates them. But it means **your MetaScope figure must go in `notebooks/`, not `data/reference/`.**
> 2. You export MetaScope figures on **your own computer** (MetaScope isn't in the Codespace). Upload the PNG back into the Codespace, drop it in `notebooks/`, and commit it — otherwise there is no record you completed the visualization.

**End-of-session commit (every session):**
```bash
git add notebooks/ log/
git commit -m "feat(moduleN): <what you did>"
git push
```

## Included Skills

Three custom skills are provided for this training:

| Skill | Description |
|-------|-------------|
| `/log` | Record a session log entry — what you did, what broke, what you learned |
| `/debug` | Structured root-cause debugging before escalating for help |
| `/explain` | Get an explanation of a concept anchored to the lab's context |

## Prerequisites

- Claude Code Pro plan account
- Basic Python (variables, loops, functions, file I/O)
- Familiarity with the terminal

## Repository Structure

```
.claude/commands/       ← /log, /debug, /explain
.devcontainer/          ← Codespace configuration
lab-context.md          ← lab context for Claude Code
notebooks/              ← all 5 training notebooks (01–05)
data/reference/         ← shared reference data and downloads
instructor/             ← rubric and instructor notes
```
