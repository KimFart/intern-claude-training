FROM mcr.microsoft.com/devcontainers/miniconda:3

# Build MEME Suite 5.5.9 from source — avoids cross-image shared library conflicts
# MPI is intentionally excluded (training use; single-node is sufficient)
RUN apt-get update && apt-get install -y \
      build-essential wget \
      libxml2-dev libxslt1-dev \
      zlib1g-dev libexpat1-dev \
      perl cpanminus \
      libhtml-template-perl libxml-simple-perl libfile-which-perl \
    && wget -q "https://meme-suite.org/meme/meme-software/5.5.9/meme-5.5.9.tar.gz" -O /tmp/meme.tar.gz \
    && tar -xzf /tmp/meme.tar.gz -C /tmp \
    && cd /tmp/meme-5.5.9 \
    && ./configure --prefix=/opt/meme --enable-opt --with-url=https://meme-suite.org 2>&1 \
    && make -j"$(nproc)" 2>&1 \
    && make install 2>&1 \
    && cd / && rm -rf /tmp/meme-5.5.9 /tmp/meme.tar.gz \
    && apt-get remove -y build-essential \
    && apt-get autoremove -y \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Configure conda channels (bioconda last added = highest priority)
RUN conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority flexible

# Create sbml env — python=3.11, all tools EXCEPT meme (compiled above)
RUN conda install -n base -c conda-forge mamba -y && \
    mamba create -n sbml -y python=3.11 \
      ipykernel jupyter notebook \
      bowtie2 samtools bedtools sra-tools \
      biopython pysam scipy pandas openpyxl matplotlib seaborn entrez-direct && \
    conda clean -afy

# Make meme + sbml env active by default — no activation needed
ENV PATH="/opt/meme/bin:/opt/conda/envs/sbml/bin:/opt/conda/bin:${PATH}"

# Install Node.js 20 and Claude Code CLI
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - && \
    apt-get install -y nodejs && \
    npm install -g @anthropic-ai/claude-code && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
