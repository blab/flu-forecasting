# Imports.
import json
import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: download_sequences, download_titers, remove_reference_strains_from_passaged_strains, summarize_model, aggregate_model_parameters, aggregate_model_accuracy, aggregate_tree_plots, aggregate_model_validation

wildcard_constraints:
    sample="(\d+|titers)",
    viruses="\d+",
    bandwidth="[0-9]*\.?[0-9]+",
    lineage="[a-z0-9]+",
    segment="[a-z]+"

# Load configuration parameters.
configfile: "config/config.json"

path_to_fauna = config["path_to_fauna"]
SEGMENTS = config["segments"]
YEAR_RANGES = config["year_ranges"]
VIRUSES = config["viruses"]
PREDICTORS = config["predictors"]
BANDWIDTHS = config["frequencies"]["bandwidths"]
MIN_LENGTH = config["min_length"]

# Construct a list of samples to build trees for including an optional tree
# where sequences with titer measurements are preferred.
NUMBER_OF_SAMPLES = config["number_of_samples"]
SAMPLES = list(range(NUMBER_OF_SAMPLES))
if config["include_titer_tree"]:
    SAMPLES += ["titers"]

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
    return expand("results/auspice/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}_tree.json", lineage="h3n2", segment=SEGMENTS, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)

include: "rules/filter_passaged_viruses.smk"
include: "rules/modular_augur_builds.smk"
include: "rules/frequency_bandwidths.smk"
include: "rules/fitness_model.smk"
include: "rules/quality_control_plots.smk"

rule all:
    input:
        "results/model_accuracy.tab",
        "results/model_parameters.tab",
        "results/model_validation.tab",
        "results/model_validation_by_bandwidth.tab",
        "results/figures/trees.pdf",
        "results/models.tab",
        "results/figures/faceted_model_fold_change.pdf",
        "results/figures/combined_model_fold_change.pdf",
        "results/figures/frequency_correlation.pdf",
        "results/figures/model_parameters.pdf",
        "results/figures/sequence_distributions.pdf",
        "results/figures/frequencies.pdf",
        _get_auspice_files,
        expand("results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/strains_metadata.tsv", lineage="h3n2", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)

rule auspice:
    input: _get_auspice_files

rule trees:
    input: rules.aggregate_tree_plots.output.trees
