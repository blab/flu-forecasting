from augur.frequency_estimators import get_pivots, timestamp_to_float
import Bio.SeqIO
import numpy as np
import pandas as pd
import pprint
import sys

from src.forecast.fitness_model import get_train_validate_timepoints

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: download_sequences, download_all_titers_by_assay, filter_metadata, filter

wildcard_constraints:
    type="(natural|simulated)",
    sample="\w+",
    timepoint="\d{4}-\d{2}-\d{2}"

# Load configuration parameters.
configfile: "config/config.json"

# Construct a list of timepoints for the requested start/end dates.
def _get_timepoints_for_build_interval(start_date, end_date, pivot_interval, min_years_per_build):
    # Find all potential timepoints.
    all_timepoints = pd.date_range(start_date, end_date, freq="%sMS" % pivot_interval)

    # Calculate date offset from the minimum years per build to find first
    # timepoint we can use to partition strains.
    offset = pd.DateOffset(years=min_years_per_build)
    first_timepoint = all_timepoints[0] + offset

    # Convert datetime instances to strings for all valid build timepoints.
    timepoints = [
        timepoint.strftime("%Y-%m-%d")
        for timepoint in all_timepoints
        if timepoint >= first_timepoint
    ]

    return timepoints

path_to_fauna = config["path_to_fauna"]

TIMEPOINT_TYPES = []
TIMEPOINT_SAMPLES = []
TIMEPOINTS = []
PREDICTOR_TYPES = []
PREDICTOR_SAMPLES = []
PREDICTORS = []
for build_type, builds_by_type in config["builds"].items():
    for sample, build in builds_by_type.items():
        # Limit Snakemake rules to active builds.
        if build["active"]:
            timepoints_for_build = _get_timepoints_for_build_interval(
                build["start_date"],
                build["end_date"],
                build["pivot_interval"],
                build["min_years_per_build"]
            )

            for timepoint in timepoints_for_build:
                TIMEPOINT_TYPES.append(build_type)
                TIMEPOINT_SAMPLES.append(sample)
                TIMEPOINTS.append(timepoint)

            for predictor in build["predictors"]:
                PREDICTOR_TYPES.append(build_type)
                PREDICTOR_SAMPLES.append(sample)
                PREDICTORS.append(predictor)

#
# Configure amino acid distance masks.
#

# Load mask configuration including which masks map to which attributes per
# lineage and segment.
masks_config = pd.read_table("config/distance_maps.tsv")

def _get_build_mask_config(wildcards):
    config = masks_config[(masks_config["lineage"] == wildcards.lineage) &
                          (masks_config["segment"] == wildcards.segment)]
    if config.shape[0] > 0:
        return config
    else:
        return None

def _get_distance_comparisons_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "compare_to"].values)

def _get_distance_attributes_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)

def _get_distance_maps_by_lineage_and_segment(wildcards):
    config = _get_build_mask_config(wildcards)
    return [
        "config/distance_maps/{wildcards.lineage}/{wildcards.segment}/{distance_map}.json".format(wildcards=wildcards, distance_map=distance_map)
        for distance_map in config.loc[:, "distance_map"].values
    ]

def _get_target_distance_earliest_date_by_wildcards(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    offset = pd.DateOffset(years=config["years_back_for_target_distance"])
    earliest_date = timepoint - offset
    return earliest_date.strftime("%Y-%m-%d")

def _get_distance_earliest_date_by_wildcards(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    season_offset = pd.DateOffset(months=config["months_for_distance_season"])
    tree_offset = pd.DateOffset(years=config["max_years_for_distances"])
    earliest_date = timepoint - season_offset - tree_offset
    return earliest_date.strftime("%Y-%m-%d")

def _get_distance_latest_date_by_wildcards(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    offset = pd.DateOffset(months=config["months_for_distance_season"])
    latest_date = timepoint - offset
    return latest_date.strftime("%Y-%m-%d")

#
# Distance functions for simulations.
#

def _get_distance_attributes_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha") &
                          (masks_config["compare_to"] != "pairwise")]
    return " ".join(config.loc[:, "attribute"].values)

def _get_pairwise_distance_attributes_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha") &
                          (masks_config["compare_to"] == "pairwise")]
    return " ".join(config.loc[:, "attribute"].values)

def _get_distance_maps_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha") &
                          (masks_config["compare_to"] != "pairwise")]
    return [
        "config/distance_maps/h3n2/ha/{distance_map}.json".format(distance_map=distance_map)
        for distance_map in config.loc[:, "distance_map"].values
    ]

def _get_pairwise_distance_maps_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha") &
                          (masks_config["compare_to"] == "pairwise")]
    return [
        "config/distance_maps/h3n2/ha/{distance_map}.json".format(distance_map=distance_map)
        for distance_map in config.loc[:, "distance_map"].values
    ]

def _get_distance_comparisons_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha") &
                          (masks_config["compare_to"] != "pairwise")]
    return " ".join(config.loc[:, "compare_to"].values)

def _get_start_date_from_range(wildcards):
    return "%s-10-01" % wildcards["year_range"].split("-")[0]

def _get_end_date_from_range(wildcards):
    return "%s-10-01" % wildcards["year_range"].split("-")[1]

def _get_predictor_list(wildcards):
    return " ".join(wildcards["predictors"].split("-"))

def _get_clock_rate_by_wildcards(wildcards):
    rates_by_lineage_and_segment = {
        ('h3n2', 'ha'): 0.0043, ('h3n2', 'na'):0.0029,
        ('h1n1pdm', 'ha'): 0.0040, ('h1n1pdm', 'na'):0.0032,
        ('vic', 'ha'): 0.0024, ('vic', 'na'):0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'):0.0013
    }

    dataset = config["datasets"][wildcards.sample]
    lineage = dataset["lineage"]
    segment = dataset["segment"]

    try:
        rate = rates_by_lineage_and_segment[(lineage, segment)]
    except KeyError:
        rate = None

    return rate

def _get_clock_rate_argument(wildcards):
    rate = _get_clock_rate_by_wildcards(wildcards)
    if rate is None:
        argument = ""
    else:
        argument = "--clock-rate %.5f" % rate

    return argument

def _get_clock_std_dev_argument(wildcards):
    rate = _get_clock_rate_by_wildcards(wildcards)
    if rate is None:
        argument = ""
    else:
        argument = "--clock-std-dev %.5f" % (0.2 * rate)

    return argument

def _get_min_date_for_augur_frequencies(wildcards):
    return timestamp_to_float(pd.to_datetime(wildcards.start))

def _get_max_date_for_augur_frequencies(wildcards):
    return timestamp_to_float(pd.to_datetime(wildcards.timepoint))

def _get_excluded_fields_arg(wildcards):
    if config.get("excluded_node_data_fields"):
        return "--excluded-fields %s" % " ".join(config["excluded_node_data_fields"])
    else:
        return ""

genes_to_translate = {
    'ha': ['SigPep', 'HA1', 'HA2'],
    'na': ['NA']
}
def gene_names(wildcards=None, segment=None):
    if segment is None:
        segment = _get_segment(wildcards)

    if segment in genes_to_translate:
        genes = genes_to_translate[segment]
    else:
        print(f"WARNING: Genes to translate are not defined for {segment}, defaulting to '{segment.upper()}'", file=sys.stderr)
        genes = [segment.upper()]

    return genes

def translations(wildcards=None, segment=None, path=None):
    genes = gene_names(wildcards, segment)
    if path is None:
        path = BUILD_TIMEPOINT_PATH

    return [path + "aa-seq_%s.fasta" % gene
            for gene in genes]

def filtered_translations(wildcards=None, segment=None, path=None):
    genes = gene_names(wildcards, segment)
    if path is None:
        path = BUILD_TIMEPOINT_PATH

    return [path + "filtered-aa-seq_%s.fasta" % gene
            for gene in genes]

#
# Define helper functions for Snakemake outputs
#

def _get_clade_model_files(wildcards):
    return expand("results/builds/{type}/{sample}/models_by_clades/{predictors}.json", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS)

def _get_distance_model_files(wildcards):
    return expand("results/builds/{type}/{sample}/models_by_distances/{predictors}.json", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS)

def _get_distance_model_errors(wildcards):
    return expand("results/builds/{type}/{sample}/annotated_models_by_distances_errors/{predictors}.tsv", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS)

def _get_distance_model_coefficients(wildcards):
    return expand("results/builds/{type}/{sample}/annotated_models_by_distances_coefficients/{predictors}.tsv", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS)

def _get_auspice_files(wildcards):
    return expand("results/auspice/flu_{type}_{sample}_{timepoint}_{filetype}.json", zip, type=TIMEPOINT_TYPES, sample=TIMEPOINT_SAMPLES, timepoint=TIMEPOINTS, filetype=["tree", "tip-frequencies"] * len(TIMEPOINTS))

#include: "rules/modular_augur_builds.smk"
#include: "rules/frequency_bandwidths.smk"
#include: "rules/fitness_model.smk"
#include: "rules/quality_control_plots.smk"
include: "rules/utils.smk"
include: "rules/datasets.smk"
include: "rules/builds.smk"

rule all:
    input:
        #_get_clade_model_files,
        _get_distance_model_files,
        _get_auspice_files,
        "results/figures/frequencies.pdf",
        "results/figures/trees.pdf"
        # "results/model_accuracy.tab",
        # "results/model_parameters.tab",
        # "results/model_validation.tab",
        # "results/model_validation_by_bandwidth.tab",
        # "results/figures/trees.pdf",
        # "results/models.tab",
        # "results/figures/faceted_model_fold_change.pdf",
        # "results/figures/combined_model_fold_change.pdf",
        # "results/figures/frequency_correlation.pdf",
        # "results/figures/model_parameters.pdf",
        # "results/figures/sequence_distributions.pdf",
        # "results/figures/frequencies.pdf",

rule clade_models:
    input: _get_clade_model_files

rule distance_models:
    input: _get_distance_model_files

rule distance_models_errors:
    input:
        errors = _get_distance_model_errors
    output:
        errors = "results/distance_model_errors.tsv"
    run:
        df = pd.concat([pd.read_csv(error_file, sep="\t") for error_file in input.errors], ignore_index=True)
        df.to_csv(output.errors, sep="\t", header=True, index=False)

rule distance_models_coefficients:
    input:
        coefficients = _get_distance_model_coefficients
    output:
        coefficients = "results/distance_model_coefficients.tsv"
    run:
        df = pd.concat([pd.read_csv(coefficient_file, sep="\t") for coefficient_file in input.coefficients], ignore_index=True)
        df.to_csv(output.coefficients, sep="\t", header=True, index=False)

rule auspice:
    input: _get_auspice_files

#rule trees:
#    input: rules.aggregate_tree_plots.output.trees
