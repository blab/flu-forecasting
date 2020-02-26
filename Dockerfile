# Start from an image with miniconda already available.
FROM continuumio/miniconda3

# Use bash as the default shell.
SHELL ["/bin/bash", "-c"]

# Create a top-level directory for all Nextstrain-related tools.
WORKDIR /nextstrain

# Install Snakemake.
RUN conda install --yes -c bioconda -c conda-forge snakemake-minimal

# Install BioPython
RUN conda install --yes -c bioconda biopython

# Install pandas
RUN conda install --yes pandas

# Download the fauna repository.
RUN git clone https://github.com/nextstrain/fauna.git

# Copy files to run the forecasting pipeline.
WORKDIR /nextstrain/flu-forecasting
COPY Snakefile .
COPY config config/
COPY data data/
COPY dist dist/
COPY envs envs/
COPY rules rules/
COPY scripts scripts/
COPY src src/
