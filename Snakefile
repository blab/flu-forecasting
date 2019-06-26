from augur.frequency_estimators import get_pivots, timestamp_to_float
import Bio.SeqIO
import numpy as np
import pandas as pd
import pprint
import sys

from src.forecast.fitness_model import get_train_validate_timepoints

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: download_sequences, download_all_titers_by_assay, filter_metadata, filter, aggregate_tree_plots

wildcard_constraints:
    sample="sample_(\d+|titers)",
    viruses="\d+",
    bandwidth="[0-9]*\.?[0-9]+",
    lineage="[a-z0-9]+",
    segment="[a-z]+[0-9]?",
    start="\d{4}-\d{2}-\d{2}",
    end="\d{4}-\d{2}-\d{2}",
    timepoint="\d{4}-\d{2}-\d{2}"

# Load configuration parameters.
configfile: "config/config.json"

path_to_fauna = config["path_to_fauna"]
LINEAGES = config["lineages"]
SEGMENTS = config["segments"]
START_DATE = config["start_date"]
END_DATE = config["end_date"]
START_DATE_TO_STANDARDIZE = config["start_date_to_standardize"]
END_DATE_TO_STANDARDIZE = config["end_date_to_standardize"]
PIVOT_INTERVAL = config["pivot_interval"]
MIN_YEARS_PER_BUILD = config["min_years_per_build"]
VIRUSES = config["viruses"]
TITER_PASSAGES = config["titers"]["passages"]
TITER_ASSAYS = config["titers"]["assays"]
PREDICTORS = config["predictors"]
BANDWIDTHS = config["frequencies"]["bandwidths"]
MIN_LENGTH = config["min_length"]

# Construct a list of samples to build trees for including an optional tree
# where sequences with titer measurements are preferred.
NUMBER_OF_SAMPLES = config["number_of_samples"]
SAMPLES = list(range(NUMBER_OF_SAMPLES))
if config["include_titer_tree"]:
    SAMPLES += ["titers"]
SAMPLES = ["sample_%s" % sample for sample in SAMPLES]

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

TIMEPOINTS = _get_timepoints_for_build_interval(START_DATE, END_DATE, PIVOT_INTERVAL, MIN_YEARS_PER_BUILD)
TRAIN_VALIDATE_TIMEPOINTS = get_train_validate_timepoints(
    TIMEPOINTS,
    config["fitness_model"]["delta_months"],
    config["fitness_model"]["training_window"]
)
#pprint.pprint(TRAIN_VALIDATE_TIMEPOINTS)
#TIMEPOINTS = TIMEPOINTS[:7]
#pprint.pprint(TIMEPOINTS)

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
                          (masks_config["segment"] == "ha")]
    return " ".join(config.loc[:, "attribute"].values)

def _get_distance_maps_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha")]
    return [
        "config/distance_maps/h3n2/ha/{distance_map}.json".format(distance_map=distance_map)
        for distance_map in config.loc[:, "distance_map"].values
    ]

def _get_distance_comparisons_for_simulations(wildcards):
    config = masks_config[(masks_config["lineage"] == "h3n2") &
                          (masks_config["segment"] == "ha")]
    return " ".join(config.loc[:, "compare_to"].values)

#
# Define helper functions.
#

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

    try:
        rate = rates_by_lineage_and_segment[(wildcards.lineage, wildcards.segment)]
    except KeyError:
        print(f"ERROR: No clock rate defined for {wildcards.lineage} and {wildcards.segment}", file=sys.stderr)
        raise

    return rate

def _get_clock_std_dev_by_wildcards(wildcards):
    return 0.2 * _get_clock_rate_by_wildcards(wildcards)

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
    if wildcards and wildcards.segment in genes_to_translate:
        genes = genes_to_translate[wildcards.segment]
    elif segment in genes_to_translate:
        genes = genes_to_translate[segment]
    else:
        print(f"WARNING: Genes to translate are not defined for {wildcards.segment}, defaulting to '{wildcards.segment.upper()}'")
        genes = [wildcards.segment.upper()]

    return genes

def translations(wildcards=None, segment=None, path=None):
    genes = gene_names(wildcards, segment)
    if path is None:
        path = BUILD_SEGMENT_PATH

    return [path + "aa-seq_%s.fasta" % gene
            for gene in genes]

#
# Define helper functions for Snakemake outputs
#

def _get_clade_model_files(wildcards):
    return expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/models_by_clades/{predictors}.json", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, predictors=PREDICTORS)

def _get_distance_model_files(wildcards):
    return expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/models_by_distances/{predictors}.json", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, predictors=PREDICTORS)

def _get_simulated_distance_model_files(wildcards):
    return expand("results/builds/simulations/{percentage}/{start}--{end}/models_by_distances/{predictors}.json", percentage=PERCENTAGE, start=START_DATE_SIMULATIONS, end=END_DATE_SIMULATIONS, predictors=PREDICTORS_SIMULATED)

def _get_auspice_files(wildcards):
    return expand("results/auspice/flu_{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}_{filetype}.json", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, timepoint=TIMEPOINTS, segment=SEGMENTS, filetype=["tree", "tip-frequencies"])

BUILD_PATH = "results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/"
BUILD_LOG_STEM = "{lineage}_{viruses}_{sample}_{start}_{end}"
BUILD_TIMEPOINT_PATH = BUILD_PATH + "timepoints/{timepoint}/"
BUILD_SEGMENT_PATH = BUILD_TIMEPOINT_PATH + "segments/{segment}/"
BUILD_SEGMENT_LOG_STEM = "{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}"

#include: "rules/filter_passaged_viruses.smk"
include: "rules/modular_augur_builds.smk"
#include: "rules/frequency_bandwidths.smk"
include: "rules/fitness_model.smk"
include: "rules/quality_control_plots.smk"
include: "rules/datasets_simulations.smk"

rule all:
    input:
        expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/tip_attributes.tsv", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE),
        expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/final_clade_frequencies.tsv", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE),
        expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/target_distances.tsv", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE),
        _get_clade_model_files,
        _get_distance_model_files,
        _get_simulated_distance_model_files,
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
        # expand("results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/strains_metadata.tsv", lineage="h3n2", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)

rule clade_models:
    input: _get_clade_model_files

rule distance_models:
    input: _get_distance_model_files

rule distance_models_simulated:
    input: _get_simulated_distance_model_files

rule auspice:
    input: _get_auspice_files

#rule trees:
#    input: rules.aggregate_tree_plots.output.trees
