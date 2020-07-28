import csv
import json
import os
import pandas as pd
from pathlib import Path
import pprint
from snakemake.logging import logger
import sys

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

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

VALIDATION_PREDICTOR_TYPES = []
VALIDATION_PREDICTOR_SAMPLES = []
VALIDATION_PREDICTOR_PREDICTORS = []

TEST_PREDICTOR_TYPES = []
TEST_PREDICTOR_SAMPLES = []
TEST_PREDICTORS = []
ACTIVE_BUILDS = config["active_builds"].split(" ")
logger.info("Active builds: " + config["active_builds"])

for build_type, builds_by_type in config["builds"].items():
    for sample, build in builds_by_type.items():
        # Limit Snakemake rules to active builds.
        if sample in ACTIVE_BUILDS:
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

            # Builds with validation samples are used for testing and not model fitting.
            if "validation_build" in build:
                validation_build = config["builds"][build_type][build["validation_build"]]
                for predictor in validation_build["predictors"]:
                    TEST_PREDICTOR_TYPES.append(build_type)
                    TEST_PREDICTOR_SAMPLES.append(sample)
                    TEST_PREDICTORS.append(predictor)
            else:
                for predictor in build["predictors"]:
                    PREDICTOR_TYPES.append(build_type)
                    PREDICTOR_SAMPLES.append(sample)
                    PREDICTORS.append(predictor)

            # Note samples for which we need validation figures.
            if "full_tree_build" in build:
                for predictor in build["validation_predictors"]:
                    VALIDATION_PREDICTOR_TYPES.append(build_type)
                    VALIDATION_PREDICTOR_SAMPLES.append(sample)
                    VALIDATION_PREDICTOR_PREDICTORS.append(predictor)

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

    sample = _get_sample_by_wildcards(wildcards)
    dataset = config["datasets"][sample]
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

def timestamp_to_float(time):
    """Convert a pandas timestamp to a floating point date.

    >>> import datetime
    >>> time = datetime.date(2010, 10, 1)
    >>> timestamp_to_float(time)
    2010.75
    >>> time = datetime.date(2011, 4, 1)
    >>> timestamp_to_float(time)
    2011.25
    >>> timestamp_to_float(datetime.date(2011, 1, 1))
    2011.0
    >>> timestamp_to_float(datetime.date(2011, 12, 1)) == (2011.0 + 11.0 / 12)
    True
    """
    return time.year + ((time.month - 1) / 12.0)

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

#
# Define helper functions for Snakemake outputs
#

def _get_clade_model_files(wildcards):
    return expand("results/builds/{type}/{sample}/models_by_clades/{predictors}.json", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS)

def _get_distance_model_files(wildcards):
    return expand("results/builds/{type}/{sample}/models_by_distances/{predictors}.json", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS) + expand("results/builds/{type}/{sample}/test_models_by_distances/{predictors}.json", zip, type=TEST_PREDICTOR_TYPES, sample=TEST_PREDICTOR_SAMPLES, predictors=TEST_PREDICTORS)

def _get_distance_model_errors(wildcards):
    return expand("results/builds/{type}/{sample}/annotated_models_by_distances_errors/{predictors}.tsv", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS) + expand("results/builds/{type}/{sample}/annotated_test_models_by_distances_errors/{predictors}.tsv", zip, type=TEST_PREDICTOR_TYPES, sample=TEST_PREDICTOR_SAMPLES, predictors=TEST_PREDICTORS)

def _get_distance_model_coefficients(wildcards):
    return expand("results/builds/{type}/{sample}/annotated_models_by_distances_coefficients/{predictors}.tsv", zip, type=PREDICTOR_TYPES, sample=PREDICTOR_SAMPLES, predictors=PREDICTORS) + expand("results/builds/{type}/{sample}/annotated_test_models_by_distances_coefficients/{predictors}.tsv", zip, type=TEST_PREDICTOR_TYPES, sample=TEST_PREDICTOR_SAMPLES, predictors=TEST_PREDICTORS)

def _get_auspice_files(wildcards):
    return expand("results/auspice/flu_{type}_{sample}_{timepoint}_{filetype}.json", zip, type=TIMEPOINT_TYPES, sample=TIMEPOINT_SAMPLES, timepoint=TIMEPOINTS, filetype=["tree", "tip-frequencies"] * len(TIMEPOINTS))

def _get_validation_figures(wildcards):
    return expand("manuscript/figures/validation_figure_{type}-{sample}-{predictors}.pdf", zip, type=VALIDATION_PREDICTOR_TYPES, sample=VALIDATION_PREDICTOR_SAMPLES, predictors=VALIDATION_PREDICTOR_PREDICTORS)

include: "rules/utils.smk"
include: "rules/datasets.smk"
include: "rules/builds.smk"

rule all:
    input:
        "results/distance_model_errors.tsv",
        "results/distance_model_coefficients.tsv",
        _get_auspice_files,
        _get_validation_figures

rule build_environment:
    conda: "envs/anaconda.python3.yaml"
    shell: "echo Environment built"

rule clade_models:
    input: _get_clade_model_files

rule distance_models:
    input: _get_distance_model_files

rule distance_models_errors:
    input:
        errors = _get_distance_model_errors
    output:
        errors = "results/distance_model_errors.tsv"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.errors} \
            --output {output.errors}
        """

rule distance_models_coefficients:
    input:
        coefficients = _get_distance_model_coefficients
    output:
        coefficients = "results/distance_model_coefficients.tsv"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/concatenate_tables.py \
            --tables {input.coefficients} \
            --output {output.coefficients}
        """

rule auspice:
    input: _get_auspice_files

rule validation_figures:
    input: _get_validation_figures

rule trees:
    input: expand(rules.aggregate_tree_plots.output.trees, zip, type=VALIDATION_PREDICTOR_TYPES, sample=VALIDATION_PREDICTOR_SAMPLES)
    output:
        trees="results/figures/trees.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

# Build figures for manuscript.

rule figure_for_model_schematic:
    input:
        tree_for_timepoint_t = "results/auspice/flu_simulated_simulated_sample_3_2029-10-01_tree.json",
        tree_for_timepoint_u = "results/auspice/flu_simulated_simulated_sample_3_2030-10-01_tree.json",
        frequencies_for_timepoint_t = "results/auspice/flu_simulated_simulated_sample_3_2029-10-01_tip-frequencies.json",
        frequencies_for_timepoint_u = "results/auspice/flu_simulated_simulated_sample_3_2030-10-01_tip-frequencies.json"
    output:
        figure = "manuscript/figures/distance-based-fitness-model.pdf"
    log:
        notebook = "logs/notebooks/plot-model-diagram.ipynb"
    conda: "envs/anaconda.python3.yaml"
    notebook:
        "notebooks/plot-model-diagram.ipynb"

rule bootstrap_analysis:
    input:
        model_distances = "results/distance_model_errors.tsv"
    output:
        output_table = "manuscript/tables/bootstrap_p_values.tsv",
        bootstrap_figure_for_simulated_sample = "manuscript/figures/bootstrap_distributions_for_simulated_sample_3.pdf",
        bootstrap_figure_for_natural_sample = "manuscript/figures/bootstrap_distributions_for_natural_sample_1_with_90_vpm_sliding.pdf",
        composite_vs_individual_model_table = "manuscript/tables/composite_vs_individual_model_comparison.tex"
    params:
        n_bootstraps = 10000
    log:
        notebook = "logs/notebooks/bootstrap-analysis.ipynb"
    conda: "envs/anaconda.python3.yaml"
    notebook:
        "notebooks/bootstrap-analysis.ipynb"

rule figures_and_tables_for_model_results:
    input:
        model_distances = "results/distance_model_errors.tsv",
        model_coefficients = "results/distance_model_coefficients.tsv",
        bootstrap_p_values = rules.bootstrap_analysis.output.output_table
    output:
        # Simulated populations
        table_for_simulated_model_selection = "manuscript/tables/simulated_model_selection.tex",
        figure_for_simulated_model_controls = "manuscript/figures/unadjusted-model-accuracy-and-coefficients-for-simulated-populations-controls.pdf",
        figure_for_simulated_individual_models = "manuscript/figures/unadjusted-model-accuracy-and-coefficients-for-simulated-populations.pdf",
        figure_for_simulated_composite_models = "manuscript/figures/unadjusted-composite-model-accuracy-and-coefficients-for-simulated-populations.pdf",

        # Natural populations
        table_for_natural_model_selection = "manuscript/tables/natural_model_selection.tex",
        table_for_natural_model_complete_selection = "manuscript/tables/complete_natural_model_selection.tex",
        figure_for_natural_epitope_vs_oracle_models = "manuscript/figures/unadjusted-composite-model-accuracy-and-coefficients-for-natural-populations-epitope-vs-oracle.pdf",
        figure_for_natural_individual_models = "manuscript/figures/unadjusted-model-accuracy-and-coefficients-for-natural-populations.pdf",
        figure_for_natural_composite_models = "manuscript/figures/best-composite-unadjusted-model-accuracy-and-coefficients-for-natural-populations.pdf",
        figure_for_natural_updated_models = "manuscript/figures/models-natural-populations-composite-with-updated-coefficients-across-test-data.pdf",

        # Cross-validation schematics
        figure_for_simulated_cross_validation = "manuscript/figures/cross-validation-for-simulated-populations.pdf",
        figure_for_natural_cross_validation = "manuscript/figures/cross-validation-for-natural-populations.pdf"
    params:
        simulated_sample = "simulated_sample_3",
        natural_sample = "natural_sample_1_with_90_vpm_sliding"
    log:
        notebook = "logs/notebooks/model-results.ipynb"
    conda: "envs/anaconda.python3.yaml"
    notebook:
        "notebooks/model-results.ipynb"

rule figure_for_vaccine_comparison:
    input:
        validation_tip_attributes = "results/builds/natural/natural_sample_1_with_90_vpm_sliding/tip_attributes_with_weighted_distances.tsv",
        test_tip_attributes = "results/builds/natural/natural_sample_1_with_90_vpm_sliding_test_tree/tip_attributes_with_weighted_distances.tsv",
        cTiter_x_ne_star_validation_forecasts_path = "results/builds/natural/natural_sample_1_with_90_vpm_sliding/forecasts_cTiter_x-ne_star.tsv",
        ne_star_lbi_validation_forecasts_path = "results/builds/natural/natural_sample_1_with_90_vpm_sliding/forecasts_ne_star-lbi.tsv",
        naive_validation_forecasts_path = "results/builds/natural/natural_sample_1_with_90_vpm_sliding/forecasts_naive.tsv",
        cTiter_x_ne_star_test_forecasts_path = "results/builds/natural/natural_sample_1_with_90_vpm_sliding_test_tree/forecasts_cTiter_x-ne_star.tsv",
        ne_star_lbi_test_forecasts_path = "results/builds/natural/natural_sample_1_with_90_vpm_sliding_test_tree/forecasts_ne_star-lbi.tsv",
        naive_test_forecasts_path = "results/builds/natural/natural_sample_1_with_90_vpm_sliding_test_tree/forecasts_naive.tsv",
        vaccines_json_path = "config/vaccines_h3n2.json"
    output:
        figure = "manuscript/figures/vaccine-comparison.pdf",
        relative_figure = "manuscript/figures/vaccine-comparison-relative-distance.pdf"
    log:
        notebook = "logs/notebooks/vaccine-strain-comparison.ipynb"
    conda: "envs/anaconda.python3.yaml"
    notebook:
        "notebooks/vaccine-strain-comparison.ipynb"

rule table_of_mutations_by_trunk_status_for_simulated_populations:
    input:
        full_tree_json = "results/auspice/flu_simulated_simulated_sample_3_full_tree_2040-10-01_tree.json",
        epitope_sites_distance_map = "config/distance_maps/h3n2/ha/luksza.json",
    output:
        table = "manuscript/tables/mutations_by_trunk_status_for_simulated_populations.tex"
    log:
        notebook = "logs/notebooks/simulated-mutations-by-trunk-status.ipynb"
    conda: "envs/anaconda.python3.yaml"
    notebook:
        "notebooks/simulated-mutations-by-trunk-status.ipynb"

rule table_of_mutations_by_trunk_status_for_natural_populations:
    input:
        full_tree_json = "results/auspice/flu_natural_natural_sample_1_with_90_vpm_sliding_full_tree_2015-10-01_tree.json",
        epitope_sites_distance_map = "config/distance_maps/h3n2/ha/luksza.json",
    output:
        table = "manuscript/tables/mutations_by_trunk_status.tex"
    log:
        notebook = "logs/notebooks/natural-mutations-by-trunk-status.ipynb"
    conda: "envs/anaconda.python3.yaml"
    notebook:
        "notebooks/natural-mutations-by-trunk-status.ipynb"

# Compile the manuscript after creating all figures and tables.
rule manuscript:
    input:
        # Tables and figures
        rules.figure_for_model_schematic.output,
        *rules.figures_and_tables_for_model_results.output,
        _get_validation_figures,
        rules.figure_for_vaccine_comparison.output,
        rules.table_of_mutations_by_trunk_status_for_simulated_populations.output,
        rules.table_of_mutations_by_trunk_status_for_natural_populations.output,

        # Manuscript text and references
        "manuscript/flu_forecasting.tex",
        "manuscript/flu_forecasting.bib",
        "manuscript/abstract.tex",
        "manuscript/main.tex",
        "manuscript/supplement.tex"
    output:
        "manuscript/flu_forecasting.pdf"
    params:
        title = "flu_forecasting"
    shell:
        """
        cd manuscript
        pdflatex -draftmode {params.title}
        bibtex {params.title}
        pdflatex -draftmode {params.title}
        pdflatex {params.title}
        """
