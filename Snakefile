# Imports.
import json
import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: clean, download_sequences, download_titers, remove_reference_strains_from_passaged_strains, summarize_model, aggregate_model_parameters, aggregate_model_accuracy, aggregate_tree_plots, aggregate_model_validation

wildcard_constraints:
    sample="(\d+|titers)",
    viruses="\d+",
    bandwidth="[0-9]*\.?[0-9]+",
    lineage="[a-z0-9]+",
    segment="[a-z]+"

# Load configuration parameters.
configfile: "config.json"

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
    return expand("auspice/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}_tree.json", lineage="h3n2", segment=SEGMENTS, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)

include: "rules/filter_passaged_viruses.smk"
include: "rules/modular_augur_builds.smk"
include: "rules/frequency_bandwidths.smk"
include: "rules/fitness_model.smk"
include: "rules/quality_control_plots.smk"

rule all:
    input:
        "model_accuracy.tab",
        "model_parameters.tab",
        "model_validation.tab",
        "model_validation_by_bandwidth.tab",
        "figures/trees.pdf",
        "models.tab",
        "figures/faceted_model_fold_change.pdf",
        "figures/combined_model_fold_change.pdf",
        "figures/frequency_correlation.pdf",
        "figures/model_parameters.pdf",
        "figures/sequence_distributions.pdf",
        "figures/frequencies.pdf",
        _get_auspice_files

rule auspice:
    input: _get_auspice_files

rule trees:
    input: "figures/trees.pdf"

rule plot_model_parameters:
    input: "model_parameters.tab"
    output: "figures/model_parameters.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/plot_parameters.py {input} {output}"

rule plot_frequency_correlation:
    input: "model_accuracy.tab"
    output:
        correlation="figures/frequency_correlation.pdf",
        mcc="figures/mcc.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/plot_accuracy.py {input} {output.correlation} {output.mcc}"

rule aggregate_model_validation:
    input: expand("model_data_frames/{year_range}/{viruses}/{predictors}/labeled_validation_{sample}.tsv", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_validation.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False, na_values="N/A") for i in input], ignore_index=True, sort=True)
        df.to_csv(output[0], sep="\t", index=False, na_rep="N/A")

rule aggregate_model_parameters:
    input: expand("model_parameters/{year_range}/{viruses}/{predictors}/{sample}.tab", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_parameters.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule aggregate_model_accuracy:
    input: expand("model_accuracy/{year_range}/{viruses}/{predictors}/{sample}.tab", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_accuracy.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False, na_rep="null")

rule aggregate_combined_model_fold_change:
    input: expand("figures/combined_model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "figures/combined_model_fold_change.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule aggregate_faceted_model_fold_change:
    input: expand("figures/faceted_model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "figures/faceted_model_fold_change.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule plot_model_fold_change:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    output:
        faceted="figures/faceted_model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf",
        combined="figures/combined_model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/plot_model_fold_change.py {input} {output.faceted} {output.combined} {wildcards.year_range} {wildcards.viruses} {wildcards.predictors} {wildcards.sample}"
