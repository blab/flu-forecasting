# Integrative prediction of seasonal influenza evolution by genotype and phenotype

## Quickstart (simulated populations only)

[Install miniconda](https://conda.io/miniconda.html) for your machine and then run the following commands.

Clone the forecasting repository.

```bash
git clone --recursive https://github.com/blab/flu-forecasting.git
cd flu-forecasting
```

Run the pipeline for simulated data.
This will first simulate influenza-like populations and then fit models to those populations.
All steps will be run locally with one CPU.

```bash
./quickstart
```

For the impatient, run locally with four CPUs.

```bash
./quickstart -j 4
```

## Run the complete analysis with simulated and natural populations

This step requires access to the Bedford lab's "fauna" database to download data for natural populations.

[Install miniconda](https://conda.io/miniconda.html) for your machine and then run the following commands.

```bash
# Clone the fauna repo
git clone https://github.com/nextstrain/fauna.git

# Clone the forecasting repo
git clone --recursive https://github.com/blab/flu-forecasting.git
cd flu-forecasting

# Create conda environment for Snakemake.
conda env create -f envs/anaconda.python3.yaml

# Set RETHINK database environment variables.
export RETHINK=

# Or run the entire pipeline on the complete input data.
./run
```

The entire pipeline is implemented with [Snakemake](https://snakemake.readthedocs.io/en/stable/), so you can run `snakemake` directly as follows to run the pipeline on your cluster.
The following example works for a SLURM-based cluster environment.
Modify the `--cluster-config` and `--drmaa` arguments to match your cluster's environment.
Change the `-j` argument to change the number of jobs to be executed simultaneously.

```bash
snakemake \
    --restart-times 3 \
    -w 60 \
    --use-conda \
    --cluster-config config/cluster-gizmo.json \
    --drmaa " -p {cluster.partition} --nodes=1 --ntasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time}" \
    --jobname "{rulename}.{jobid}.sh"
    -j 20
```

By default, this pipeline will download HA sequences and public titers for A/H3N2 with fauna, run augur prepare and process, run the fitness model for all defined combinations of predictors, and generate the output tables and figures described below.

## Configuration

Model builds are parameterized by the contents of `config.json`.
The following parameters are required to specify builds.

| Parameter | Description | Example |
|:---:|:---:|:---:|
| year_ranges | list of year intervals to build trees for | `["2000-2015", "2005-2017"]` |
| viruses | scalar or list of number of viruses to sample per month | `[20, 92]` |
| number_of_samples | number of trees to build for each combination of year ranges and virus sampling densities | `5` |
| predictors | list of fitness predictors to fit a model for; multiple predictors specified as hyphen-delimited lists | `["null", "ep-cTiterSub-dms"]` |

Trees will be built for all combination of year ranges, virus sampling densities, and number of samples.
Models will be built for all combination of predictors and trees.

## Inputs

The only inputs currently are the configuration file, `config.json`, and the data downloaded from fauna.

## Outputs

| Filename | Contents |
|:---:|:---:|
| models.tab | raw data table of initial, observed, and predicted frequencies by clade, timepoint, and model |
| model_accuracy.tab | summary table of growth correlation and clade error for all models |
| model_parameters.tab | summary table of beta and std dev params for each predictor in each model with one or more predictors |
| trees.pdf | BALTIC-style visualization of trees built by augur's process step |
| model_fold_change.pdf | Correlation of observed and predicted growth rate for all models |
