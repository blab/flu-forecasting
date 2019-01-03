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

rule aggregate_frequency_plots:
    input: expand("figures/frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "figures/frequencies.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule aggregate_sequence_distribution_plots:
    input: expand("figures/sequence_distributions/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "figures/sequence_distributions.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule aggregate_tree_plots:
    input: expand("figures/trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.pdf", segment=SEGMENTS, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "figures/trees.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

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

rule aggregate_models:
    input: expand("models/{year_range}/{viruses}/{predictors}/{sample}.tab", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "models.tab"
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

rule label_model_validation_table:
    input: "model_data_frames/{year_range}/{viruses}/{predictors}/validation_{sample}.tsv"
    output: "model_data_frames/{year_range}/{viruses}/{predictors}/labeled_validation_{sample}.tsv"
    run:
        df = pd.read_table(input[0])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample
        df.to_csv(output[0], sep="\t", header=True, index=False, na_rep="N/A")

rule convert_model_json_to_tsv:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    output: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    run:
        with open(input[0], "r") as fh:
            model_json = json.load(fh)

        df = pd.DataFrame(model_json["test_data"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample

        df.to_csv(output[0], sep="\t", header=True, index=False)

rule summarize_model:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    output:
        accuracy="model_accuracy/{year_range}/{viruses}/{predictors}/{sample}.tab",
        parameters="model_parameters/{year_range}/{viruses}/{predictors}/{sample}.tab"
    run:
        with open(input[0], "r") as fh:
            model = json.load(fh)

        accuracy = model["accuracy"]
        df = pd.DataFrame({
            "year_range": [wildcards.year_range],
            "viruses": [wildcards.viruses],
            "predictors": [wildcards.predictors],
            "sample": [wildcards.sample],
            "correlation_rel": [accuracy["correlation_rel"]],
            "mcc": [accuracy["mcc"]],
            "clade_error": [accuracy["clade_error"]]
        })
        df.to_csv(output["accuracy"], sep="\t", index=False, na_rep="NaN")

        df = pd.DataFrame(model["params"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample
        df.to_csv(output["parameters"], sep="\t", index=False)

rule run_fitness_model:
    input:
        ha_tree="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        ha_metadata="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_meta.json",
        ha_sequences="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_seq.json",
        #na_tree="auspice/flu_h3n2_na_{year_range}y_{viruses}v_{sample}_tree.json",
        frequencies="frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        titers=rules.download_titers.output.titers,
        dms="data/dms-h3n2-preferences-rescaled.csv"
    output:
        model="models/{year_range}/{viruses}/{predictors}/{sample}.json",
        tip_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/tips_{sample}.tsv",
        clade_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/clades_{sample}.tsv",
        validation_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/validation_{sample}.tsv"
    params:
        predictor_list=_get_predictor_list,
        min_freq=config["fitness_model"]["min_freq"],
        max_freq=config["fitness_model"]["max_freq"]
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.txt"
    log: "logs/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.log"
    shell: "python fit_model.py {input.ha_tree} {input.ha_metadata} {input.ha_sequences} {input.frequencies} {output.model} {params.predictor_list} --titers {input.titers} --dms {SNAKEMAKE_DIR}/{input.dms} --tip-data-frame {output.tip_data_frame} --clade-data-frame {output.clade_data_frame} --validation-data-frame {output.validation_data_frame} --min-freq {params.min_freq} --max-freq {params.max_freq} -v &> {log}"

rule plot_frequencies:
    input: "frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    output: "figures/frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/plot_frequency_trajectories.py {input} {output}"

rule estimate_frequencies:
    message:
        """
        Estimating frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        weights="data/region_weights.json"
    output: "frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=config["frequencies"]["pivot_frequency"],
        start_date=_get_start_date_from_range,
        end_date=_get_end_date_from_range
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/estimate_frequencies_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/estimate_frequencies_{year_range}y_{viruses}v_{sample}.log"
    shell: """python scripts/frequencies.py {input.tree} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--pivot-frequency {params.pivot_frequency} \
--start-date {params.start_date} \
--end-date {params.end_date} \
--weights {input.weights} \
--weights-attribute region \
--include-internal-nodes \
--censored &> {log}"""

rule plot_sequences_by_date:
    input: "auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json"
    output: "figures/sequence_distributions/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/plot_sequence_distribution.py {input} {output}"

rule plot_tree:
    input: "auspice/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.json",
    output: "figures/trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.pdf"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/plot_tree_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/plot_tree_{segment}_{year_range}y_{viruses}v_{sample}.log"
    shell: """cd dist/augur/scripts && python plot_tree.py {SNAKEMAKE_DIR}/{input} {SNAKEMAKE_DIR}/{output} &> {SNAKEMAKE_DIR}/{log}"""
