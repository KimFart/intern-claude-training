#!/bin/bash
set -e

echo "==> Setting up SBML intern training environment..."

# Source conda for this script only (does not affect terminals)
source /opt/conda/etc/profile.d/conda.sh

# Configure channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install mamba for faster resolution
echo "==> Installing mamba..."
conda install -n base -c conda-forge mamba -y

# Create sbml env with all packages at once
echo "==> Creating sbml conda environment..."
mamba create -n sbml -y python=3.11 \
  ipykernel jupyter notebook \
  bowtie2 samtools bedtools sra-tools meme \
  biopython pysam scipy pandas matplotlib seaborn entrez-direct

# Register kernel using direct path — no conda run, no conda activate
echo "==> Registering Jupyter kernel..."
/opt/conda/envs/sbml/bin/python -m ipykernel install --user \
  --name sbml --display-name "Python (sbml)"

# Install Claude Code CLI
echo "==> Installing Claude Code CLI..."
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs
npm install -g @anthropic-ai/claude-code

echo "==> Setup complete."
