import pandas as pd
import pprint

from src.forecast.fitness_model import get_train_validate_timepoints

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: download_sequences, download_titers, remove_reference_strains_from_passaged_strains, summarize_model, aggregate_model_parameters, aggregate_model_accuracy, aggregate_tree_plots, aggregate_model_validation

wildcard_constraints:
    sample="sample_(\d+|titers)",
    viruses="\d+",
    bandwidth="[0-9]*\.?[0-9]+",
    lineage="[a-z0-9]+",
    segment="[a-z]+",
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
TRAIN_VALIDATE_TIMEPOINTS = get_train_validate_timepoints(TIMEPOINTS, 1, 1)
pprint.pprint(TRAIN_VALIDATE_TIMEPOINTS)
TIMEPOINTS = TIMEPOINTS[:3]

#
# Configure amino acid distance masks.
#

# Load mask configuration including which masks map to which attributes per
# lineage and segment.
masks_config = pd.read_table("config/mask_config.tsv")

def _get_build_mask_config(wildcards):
    config = masks_config[(masks_config["lineage"] == wildcards.lineage) &
                          (masks_config["segment"] == wildcards.segment)]
    if config.shape[0] > 0:
        return config
    else:
        return None

def _get_mask_attribute_names_by_wildcards(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)

def _get_mask_names_by_wildcards(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "mask"].values)

#
# Define helper functions.
#

def _get_start_date_from_range(wildcards):
    return "%s-10-01" % wildcards["year_range"].split("-")[0]

def _get_end_date_from_range(wildcards):
    return "%s-10-01" % wildcards["year_range"].split("-")[1]

def _get_predictor_list(wildcards):
    return " ".join(wildcards["predictors"].split("-"))

def _get_distance_attribute_names_by_segment(wildcards):
    return " ".join(config["distance_parameters_by_segment"][wildcards.segment]["attributes"])

def _get_distance_mask_names_by_segment(wildcards):
    return " ".join(config["distance_parameters_by_segment"][wildcards.segment]["masks"])

def _get_clock_rate_by_wildcards(wildcards):
    rate = {
        ('h3n2', 'ha'): 0.0043, ('h3n2', 'na'):0.0029,
        ('h1n1pdm', 'ha'): 0.0040, ('h1n1pdm', 'na'):0.0032,
        ('vic', 'ha'): 0.0024, ('vic', 'na'):0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'):0.0013
    }
    return rate[(wildcards.lineage, wildcards.segment)]

def _get_auspice_files(wildcards):
    return expand("results/auspice/flu_{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}_tree.json", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, timepoint=TIMEPOINTS, segment=SEGMENTS)

BUILD_PATH = "results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/"
BUILD_TIMEPOINT_PATH = BUILD_PATH + "timepoints/{timepoint}/"
BUILD_SEGMENT_PATH = BUILD_TIMEPOINT_PATH + "segments/{segment}/"
BUILD_SEGMENT_LOG_STEM = "{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}"

#include: "rules/filter_passaged_viruses.smk"
include: "rules/modular_augur_builds.smk"
#include: "rules/frequency_bandwidths.smk"
include: "rules/fitness_model.smk"
include: "rules/quality_control_plots.smk"

rule all:
    input:
        expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/tip_attributes.tsv", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE),
        expand("results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/tips_to_clades.tsv", lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE),
        _get_auspice_files,
        "results/figures/frequencies.pdf"
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

rule auspice:
   input: _get_auspice_files

#rule trees:
#    input: rules.aggregate_tree_plots.output.trees
