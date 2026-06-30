#!/bin/bash
set -e

echo "==> Setting up SBML intern training environment..."

# Make conda commands available in this non-interactive script
source /opt/conda/etc/profile.d/conda.sh

# Configure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install mamba for faster package resolution
echo "==> Installing mamba..."
conda install -n base -c conda-forge mamba -y

# Phase 1: Create env with Python + Jupyter immediately so VS Code has a valid interpreter
echo "==> Phase 1: Creating sbml env with Python + Jupyter..."
mamba create -n sbml -y python=3.11 ipykernel jupyter notebook

# Register the kernel right away so VS Code can use it
conda run -n sbml python -m ipykernel install --user --name sbml --display-name "Python (sbml)"

echo "==> Python (sbml) kernel ready. VS Code can now connect."

# Phase 2: Install bioinformatics tools
echo "==> Phase 2: Installing bioinformatics tools (this takes a while)..."
mamba install -n sbml -y \
  bowtie2 \
  samtools \
  bedtools \
  sra-tools \
  meme \
  biopython \
  pysam \
  scipy \
  pandas \
  matplotlib \
  seaborn \
  entrez-direct

# Auto-activate sbml environment in new shells
echo "conda activate sbml" >> ~/.bashrc

# Install Claude Code CLI (requires Node.js)
echo "==> Installing Claude Code CLI..."
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs
npm install -g @anthropic-ai/claude-code

echo ""
echo "==> Environment fully ready. Run 'claude' in the terminal and authenticate with your Pro plan account."
