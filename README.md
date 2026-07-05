# SBML Intern Claude Code Training

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/sbml-lab/intern-claude-training)

> The Codespace will install all required tools automatically (bowtie2, samtools, MEME, Biopython, Claude Code CLI). First launch takes ~5 minutes.
>
> You will need to authenticate Claude Code with a **Pro plan** account after launch: run `claude` in the terminal.

## Overview

A 6-week training curriculum for SBML Lab (KAIST GSEB) interns. Covers Python, GFF/pandas data processing, NGS alignment, motif analysis, and an independent mini-project — all using Claude Code as the primary tool.

## Curriculum

| Module | Topic |
|--------|-------|
| 1 | What is Claude Code |
| 2 | GFF parsing + pandas with Claude Code |
| 3 | NGS alignment pipeline with Claude Code |
| 4 | Papers, MEME, and Biopython with Claude Code |
| 5 (Weeks 5–6) | Independent mini-project |

## Standing Rules

All interns follow these rules throughout the training:

1. `/log` at the end of every session
2. `/debug` before asking for help when something breaks
3. `/explain [concept]` before Googling
4. Plan mode (Shift+Tab twice) before writing any pipeline

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
