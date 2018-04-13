# Seasonal influenza fitness predictions

## Quick start

[Install miniconda](https://conda.io/miniconda.html) for your machine and then run the following commands.

```bash
# Create conda environment for Snakemake.
conda env create -f envs/anaconda.python3.yaml

# Load Snakemake environment.
source activate janus_python3

# Run pipeline.
# The first run will take some time while the conda environment is created.
snakemake --use-conda
```
