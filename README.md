# Integrative prediction of seasonal influenza evolution by genotype and phenotype

## Quick start

[Install miniconda](https://conda.io/miniconda.html) for your machine and then run the following commands.

```bash
# Clone the repo
git clone --recursive https://github.com/blab/flu-forecasting.git
cd flu-forecasting

# Create conda environment for Snakemake.
conda env create -f envs/anaconda.python3.yaml

# Set RETHINK database environment variables.
export RETHINK=

# Run pipeline on a minimal dataset to confirm everything works.
# The first run will take some time while Snakemake creates the conda environments it needs.
./quickstart

# Or run the entire pipeline on the complete input data.
./run
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
