# SBML Lab — Context for Claude Code

This file is the lab context for Claude Code skills and prompts.
When using `/explain`, `/debug`, or any Claude Code interaction, this is the background Claude should assume.

---

## Lab Identity

- **Lab**: SBML (Systems Biology & Metabolic Engineering Lab)
- **PI**: Prof. Donghyuk Kim (김동혁)
- **Institution**: KAIST GSEB (Graduate School of Engineering Biology)
- **Website**: sbml-lab.ai

---

## Organisms

| Organism | Strain | Use |
|---|---|---|
| *E. coli* K-12 | MG1655 | Primary model organism |

Reference genome for *E. coli* K-12 MG1655: **NC_000913** (NCBI)

---

## Experimental Techniques

| Technique | Purpose |
|---|---|
| ChIP-exo / ChIP-mini | Protein-DNA binding site mapping (transcription factors) |
| RNA-seq | Differential gene expression |

---

## Bioinformatics Stack

### Alignment
- **Aligner**: bowtie2
- **Reference indexing**: `bowtie2-build`
- **Typical single-end command**: `bowtie2 -x <index> -U <reads.fastq> -S <out.sam>`
- **Typical paired-end command**: `bowtie2 -x <index> -1 <R1.fastq> -2 <R2.fastq> -S <out.sam> -X 1000`

### Post-alignment
- **SAM → BAM**: `samtools view -bS`
- **Sort**: `samtools sort`
- **Index**: `samtools index`
- **Coverage/GFF**: `makegff.py` (lab custom script — converts BAM to GFF format for MetaScope)

### Data Retrieval
- Source: NCBI GEO / SRA
- Download: `fastq-dump --split-files <SRR_accession>`

### File Formats
| Format | Description |
|---|---|
| FASTQ | Raw sequencing reads |
| SAM | Alignment (text) |
| BAM | Alignment (binary, sorted + indexed for downstream use) |
| GFF | Genomic feature format — used as input for MetaScope |
| FASTA | Reference genome / sequences |

### GFF Format (9 columns, tab-delimited)
```
seqname  source  feature  start  end  score  strand  frame  attribute
```
- 1-based coordinates
- Strand: `+` or `-`
- The annotation GFF uses `gene` as the feature type; `makegff.py` output (aligned reads) uses `read` as the feature type

### Python Libraries
- `pandas` — DataFrame operations on GFF / count tables
- `Biopython` (`Bio.SeqIO`, `Bio.Seq`) — FASTA parsing, sequence extraction
- `matplotlib`, `seaborn` — visualization

### Motif Analysis
- Tool: MEME Suite (`meme` command)
- Input: multi-FASTA of binding site sequences (typically ±50bp around ChIP peak summits)
- Key parameters: `-dna -mod zoops -nmotifs 5 -minw 6 -maxw 20`
- Output: `meme.html` (visual), `meme.txt` (parseable)
- Interpretation: E-value < 0.05 is typically considered significant

---

## Visualization

- **MetaScope**: the lab's own genome browser (cross-platform, Avalonia)
- Input: GFF files produced by `makegff.py`
- Used to visually inspect ChIP-exo peaks, RNA-seq coverage, and multi-dataset overlays

---

## Naming Conventions

- GFF output files: `<sample>_<condition>.gff`
- BAM files: `<sample>_sorted.bam`
- Index files: `<sample>_sorted.bam.bai`
- FASTQ from SRA: `<SRR_accession>_1.fastq`, `<SRR_accession>_2.fastq` (paired-end)

---

## Common Workflows

### ChIP-exo pipeline (single-end)
```
fastq-dump <SRR> → bowtie2 align → samtools sort/index → makegff.py → MetaScope
```

### RNA-seq pipeline (paired-end)
```
fastq-dump --split-files <SRR> → bowtie2 paired-end → samtools sort/index → makegff.py → MetaScope
```

### Motif discovery
```
Extract peak sequences (±50bp) with Biopython → write multi-FASTA → run MEME → interpret motif
```
