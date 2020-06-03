# Integrative prediction of seasonal influenza evolution by genotype and phenotype

## Installation

The software for this pipeline is available from GitHub + conda or a Docker image.

### Installation with GitHub and conda

[Install miniconda](https://conda.io/miniconda.html).
Clone the forecasting repository.

```bash
git clone https://github.com/blab/flu-forecasting.git
cd flu-forecasting
```

Create and activate a conda environment for the pipeline.

```bash
conda env create -f envs/anaconda.python3.yaml
conda activate flu_forecasting
```

### Installation with Docker

Install [Docker](https://www.docker.com/).
Download the image for the forecasting pipeline.

```
docker pull nextstrain/flu-forecasting:latest
```

Start a shell inside a container from this image.

```
docker run -i -t flu-forecasting:latest
```

## Quickstart (sparse simulated populations only)

Run the pipeline for sparse simulated data.
This will first simulate influenza-like populations and then fit models to those populations.
Inspect all steps to be executed by the pipeline with a dryrun.

```bash
snakemake --dryrun --config active_builds='simulated_sample_1'
```

Run the pipeline locally with four jobs (or cores) at once.

```bash
snakemake --config active_builds='simulated_sample_1' -j 4
```

Always specify a value for `-j`, to limit the number of cores available to the simulator.
If no limit is provided, the Java-based simulator will attempt to use all available cores and may cause headaches for you or your cluster's system administrator.

## Run the complete analysis with simulated and natural populations

This step requires access to the Bedford lab's "fauna" database to download data for natural populations.
Clone the fauna repository in the parent directory of the `flu-forecasting` directory.

```bash
git clone https://github.com/nextstrain/fauna.git ../fauna
```

Set environment variables to connect to the database.

```bash
export RETHINK_HOST=yourhost
export RETHINK_AUTH_KEY=yourkey
```

Run the entire pipeline locally with four simultaneous jobs.

```bash
snakemake -j 4
```

Alternately, follow [Snakemake documentation to distribute the entire pipeline to your cloud or cluster accounts](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html).
The following is an example of how to distribute the pipeline on a SLURM-based cluster using a Snakemake profile.

```bash
snakemake --profile profiles/slurm-drmaa
```

You can also run just one of the natural builds as follows, to confirm your environment is configured properly.

```bash
snakemake --config active_builds='natural_sample_1_with_90_vpm_sliding' -j 4
```

## Configuration

Analyses are parameterized by the contents of `config.json`.
Models are fit to annotated data frames created for one or more "builds" from one or more "datasets".
Datasets and builds are decoupled to allow multiple builds from a single dataset.
Builds are split into "simulated" and "natural" such that each entry in one of these categories is a dictionary of build settings indexed by a build name.
The list of active builds is determined by the space-delimited values in the `active_builds` top-level key of the configuration.

## Inputs

The only inputs currently are the configuration file, `config.json`, and the data downloaded from fauna.

## Outputs

| Filename | Contents |
|:---:|:---:|
| results/distance_model_errors.tsv | table of model errors per build, timepoint, and model |
| results/distance_model_coefficients.tsv | table of model coefficients per build, timepoint, and model |
