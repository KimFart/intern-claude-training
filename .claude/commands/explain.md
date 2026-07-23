---
name: explain
description: Explains any tool, concept, or file format anchored to the SBML lab's context. Better than Googling for lab-specific questions.
---

# Explain (Lab-Anchored)

## Purpose

Give a clear, practical explanation of a bioinformatics concept, tool, or file format — grounded in what the SBML lab actually uses.

## Lab Context

Always assume this environment when explaining:
- Organism: *E. coli* K-12 MG1655
- Aligner: bowtie2
- Post-alignment: samtools (sort, index, view)
- Coverage: makegff.py → GFF → MetaScope
- Motif discovery: MEME Suite
- Python: pandas, Biopython
- Experiment types: ChIP-exo (single-end), RNA-seq (paired-end)

## Process

1. Read the user's concept or question.

2. Explain it in three parts:
   - **What it is**: one sentence definition
   - **How it works**: 2–4 sentences on the mechanism or logic
   - **In our lab**: how this concept applies specifically to the lab's tools, files, or workflows

3. If relevant, include a concrete example using the lab's actual commands, file formats, or organism.

4. Keep the total explanation under 200 words. If the topic needs more depth, ask if the user wants to go deeper.

## Example Output Format

**[Concept]**

*What it is*: ...

*How it works*: ...

*In our lab*: ...

## Notes

- Do not give generic textbook explanations. Anchor every explanation to the lab's actual stack.
- If the concept doesn't directly relate to the lab's work, still use the lab's tools as analogies where helpful.
- If the user extends this skill with additional details from `lab-context.md`, incorporate them into the "In our lab" section.
