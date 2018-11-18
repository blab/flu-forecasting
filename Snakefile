# Imports.
import json
import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: clean, download_sequences_and_titers, remove_reference_strains_from_passaged_strains, augur_prepare, summarize_model, aggregate_model_parameters, aggregate_model_accuracy, aggregate_tree_plots

ruleorder: augur_prepare_with_titers > augur_prepare

wildcard_constraints:
    sample="(\d+|titers)",
    viruses="\d+"

# Load configuration parameters.
configfile: "config.json"

SEGMENTS = config["segments"]
YEAR_RANGES = config["year_ranges"]
VIRUSES = config["viruses"]
PREDICTORS = config["predictors"]

# Construct a list of samples to build trees for including an optional tree
# where sequences with titer measurements are preferred.
NUMBER_OF_SAMPLES = config["number_of_samples"]
SAMPLES = list(range(NUMBER_OF_SAMPLES))
if config["include_titer_tree"]:
    SAMPLES += ["titers"]

def _get_start_date_from_range(wildcards):
    return "%s-10-01" % wildcards["year_range"].split("-")[0]

def _get_end_date_from_range(wildcards):
    return "%s-04-01" % wildcards["year_range"].split("-")[1]

def _get_predictor_list(wildcards):
    return " ".join(wildcards["predictors"].split("-"))

def _get_distance_attribute_names_by_segment(wildcards):
    return " ".join(config["distance_parameters_by_segment"][wildcards.segment]["attributes"])

def _get_distance_mask_names_by_segment(wildcards):
    return " ".join(config["distance_parameters_by_segment"][wildcards.segment]["masks"])

rule all:
    input:
        "model_accuracy.tab",
        "model_parameters.tab",
        "figures/trees.pdf",
        "models.tab",
        "figures/faceted_model_fold_change.pdf",
        "figures/combined_model_fold_change.pdf",
        "figures/frequency_correlation.pdf",
        "figures/model_parameters.pdf",
        "figures/sequence_distributions.pdf",
        "figures/frequencies.pdf"

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

rule convert_model_json_to_tsv:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    output: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    run:
        with open(input[0], "r") as fh:
            model_json = json.load(fh)

        df = pd.DataFrame(model_json["data"])
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
        ha_tree="trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        ha_metadata="metadata/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_meta.json",
        ha_sequences="sequences/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_seq.json",
        #na_tree="trees/flu_h3n2_na_{year_range}y_{viruses}v_{sample}_tree.json",
        frequencies="frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        titers="dist/fauna/data/h3n2_who_hi_cell_titers.tsv",
        dms="data/dms-h3n2-preferences-rescaled.csv"
    output:
        model="models/{year_range}/{viruses}/{predictors}/{sample}.json",
        tip_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/tips_{sample}.tsv",
        clade_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/clades_{sample}.tsv"
    params:
        predictor_list=_get_predictor_list,
        min_freq=config["fitness_model"]["min_freq"],
        max_freq=config["fitness_model"]["max_freq"]
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.txt"
    log: "logs/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.log"
    shell: "python fit_model.py {input.ha_tree} {input.ha_metadata} {input.ha_sequences} {input.frequencies} {output.model} {params.predictor_list} --titers {input.titers} --dms {SNAKEMAKE_DIR}/{input.dms} --tip-data-frame {output.tip_data_frame} --clade-data-frame {output.clade_data_frame} --min-freq {params.min_freq} --max-freq {params.max_freq} -v &> {log}"

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
        tree="trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
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
--include-internal-nodes &> {log}"""

rule plot_sequences_by_date:
    input: "trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json"
    output: "figures/sequence_distributions/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/plot_sequence_distribution.py {input} {output}"

rule plot_tree:
    input: "trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.json",
    output: "figures/trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.pdf"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/plot_tree_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/plot_tree_{segment}_{year_range}y_{viruses}v_{sample}.log"
    shell: """cd dist/augur/scripts && python plot_tree.py {SNAKEMAKE_DIR}/{input} {SNAKEMAKE_DIR}/{output} &> {SNAKEMAKE_DIR}/{log}"""

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        metadata = "builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.tsv",
        branch_lengths = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/branch_lengths.json",
        traits = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/traits.json",
        nt_muts = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/nt_muts.json",
        aa_muts = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aa_muts.json",
        distances = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/distances.json",
        lbi="builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/lbi.json",
        clades = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/clades.json",
        colors = "dist/augur/builds/flu/colors.tsv",
        auspice_config = "config/auspice_config.json"
    output:
        auspice_tree = "trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.json",
        auspice_metadata = "metadata/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_meta.json",
        auspice_sequence = "sequences/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_seq.json",
    params:
        geography_traits = "region",
        panels = "tree entropy"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.distances} \
                        {input.lbi} \
                        {input.clades} \
            --colors {input.colors} \
            --geography-traits {params.geography_traits} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_metadata} \
            --output-sequence {output.auspice_sequence} \
            --panels {params.panels}
        """

rule lbi:
    message: "Calculating LBI with tau={params.tau} and window={params.window}"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        branch_lengths = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/branch_lengths.json",
    output:
        lbi = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/lbi.json"
    conda: "envs/anaconda.python2.yaml"
    params:
        tau = 0.5,
        window = 0.1
    shell:
        """
        python scripts/lbi.py \
            {input.tree} \
            {input.branch_lengths} \
            {output.lbi} \
            --attribute-names lbi \
            --tau {params.tau} \
            --window {params.window}
        """

rule clades:
    message: "Annotating clades"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        nt_muts = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/nt_muts.json",
        aa_muts = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aa_muts.json",
        clade_definitions = "clades.tsv"
    output:
        "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/clades.json"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nt_muts} {input.aa_muts} \
            --clades {input.clade_definitions} \
            --output {output}
        """

rule distances:
    message: "Calculating amino acid distances"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        sequences = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/sequences.json",
        masks = "dist/augur/builds/flu/metadata/{segment}_masks.tsv"
    output:
        distances = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/distances.json"
    params:
        attribute_names = _get_distance_attribute_names_by_segment,
        mask_names = _get_distance_mask_names_by_segment
    conda: "envs/anaconda.python2.yaml"
    shell:
        """
        python scripts/distance.py \
            {input.tree} \
            {input.sequences} \
            {output.distances} \
            --masks {input.masks} \
            --attribute-names {params.attribute_names} \
            --mask-names {params.mask_names}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        metadata = "builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.tsv"
    output:
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/traits.json",
    params:
        columns = "region country"
    conda: "envs/anaconda.python3.yaml"
    benchmark: "benchmarks/traits_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule reconstruct_sequences:
    message: "Reconstructing sequences from tree and root sequences"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        nt_muts="builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/nt_muts.json",
        aa_muts="builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aa_muts.json"
    output:
        "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/sequences.json"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/reconstruct_sequences.py {input.tree} {input.nt_muts} {input.aa_muts} {output}"

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/nt_muts.json",
        reference = "dist/augur/builds/flu/metadata/h3n2_{segment}_outgroup.gb"
    output:
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aa_muts.json"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        alignment = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aligned.fasta"
    output:
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/nt_muts.json"
    params:
        inference = "joint"
    conda: "envs/anaconda.python3.yaml"
    benchmark: "benchmarks/ancestral_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree_raw.nwk",
        alignment = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aligned.fasta",
        metadata = "builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.tsv"
    output:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/branch_lengths.json"
    params:
        coalescent = 0.03,
        clock_filter_iqd = 4
    conda: "envs/anaconda.python3.yaml"
    benchmark: "benchmarks/refine_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aligned.fasta"
    output:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree_raw.nwk"
    conda: "envs/anaconda.python3.yaml"
    shadow: "minimal"
    benchmark: "benchmarks/tree_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --nthreads {threads}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = "builds/data/filtered/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.fasta",
        reference = "dist/augur/builds/flu/metadata/h3n2_{segment}_outgroup.gb"
    output:
        alignment = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aligned.fasta"
    conda: "envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    threads: 4
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps \
            --nthreads {threads}
        """

rule filter:
    message:
        """
        Filtering sequences to exclude H1N1 and swine samples
        """
    input:
        sequences = "builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.fasta",
        metadata = "builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.tsv",
        excluded = "outliers/filtered_strains.txt"
    output:
        sequences = "builds/data/filtered/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.fasta"
    conda: "envs/anaconda.python3.yaml"
    benchmark: "benchmarks/filter_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.excluded} \
            --output {output.sequences}
        """

rule convert_prepared_json_to_metadata_and_sequences:
    input: "dist/augur/builds/flu/prepared/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.json"
    output:
        sequences="builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.fasta",
        metadata="builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.tsv"
    conda: "envs/anaconda.python2.yaml"
    log: "logs/convert_prepared_json_{segment}_{year_range}y_{viruses}v_{sample}.log"
    shell: "cd dist/augur/scripts && python prepared_json_to_fasta.py --metadata {SNAKEMAKE_DIR}/{output.metadata} {SNAKEMAKE_DIR}/{input} > {SNAKEMAKE_DIR}/{output.sequences} 2> {SNAKEMAKE_DIR}/{log}"

rule augur_prepare_with_titers:
    input:
        ha_sequences="dist/fauna/data/h3n2_ha_unpassaged.fasta",
        na_sequences="dist/fauna/data/h3n2_na.fasta",
        titers="dist/fauna/data/h3n2_who_hi_cell_titers.tsv"
    output:
        "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_titers.json",
        "dist/augur/builds/flu/prepared/flu_h3n2_na_{year_range}y_{viruses}v_titers.json"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/augur_prepare_{year_range}y_{viruses}v_titers.txt"
    log: "logs/augur_prepare_{year_range}y_{viruses}v_titers.log"
    params: start_date=_get_start_date_from_range, end_date=_get_end_date_from_range
    shell: """cd dist/augur/builds/flu && python flu.prepare.py -v {wildcards.viruses} \
  --sequences ../../../../{input.ha_sequences} ../../../../{input.na_sequences} \
  --file_prefix flu_h3n2_*segment*_{wildcards.year_range}y_{wildcards.viruses}v_titers --lineage h3n2 --segments ha na --time_interval {params.start_date} {params.end_date} -r 12y --sampling even --titers {SNAKEMAKE_DIR}/{input.titers} &> {SNAKEMAKE_DIR}/{log}"""

rule augur_prepare:
    input:
        ha_sequences="dist/fauna/data/h3n2_ha_unpassaged.fasta",
        na_sequences="dist/fauna/data/h3n2_na.fasta"
    output:
        "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        "dist/augur/builds/flu/prepared/flu_h3n2_na_{year_range}y_{viruses}v_{sample}.json"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/augur_prepare_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/augur_prepare_{year_range}y_{viruses}v_{sample}.log"
    params: start_date=_get_start_date_from_range, end_date=_get_end_date_from_range
    shell: """cd dist/augur/builds/flu && python flu.prepare.py -v {wildcards.viruses} \
  --sequences ../../../../{input.ha_sequences} ../../../../{input.na_sequences} \
  --file_prefix flu_h3n2_*segment*_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample} --lineage h3n2 --segments ha na --time_interval {params.start_date} {params.end_date} -r 12y --sampling even &> {SNAKEMAKE_DIR}/{log}"""

rule filter_passaged_viruses:
    message:
        """
        Filtering to exclude passaged viruses
        """
    input:
        sequences = "dist/fauna/data/h3n2_ha.fasta",
        metadata = "filtering_metadata.tsv",
        excluded = "outliers/non_reference_passaged_strains.txt"
    output:
        sequences = "dist/fauna/data/h3n2_ha_unpassaged.fasta"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.excluded} \
            --output {output.sequences}
        """

rule build_simple_metadata_for_filtering:
    input: "dist/fauna/data/h3n2_ha.fasta"
    output: "filtering_metadata.tsv"
    shell: """grep "^>" {input} | sed 's/>//' | sort | uniq | awk '{{ if (NR == 1) {{ print "strain" }} print }}' > {output}"""

rule remove_reference_strains_from_passaged_strains:
    message: "Removing reference strains from list of passaged strains"
    input:
        passaged="outliers/passaged_strains.txt",
        references="data/reference_viruses.txt"
    output: "outliers/non_reference_passaged_strains.txt"
    run:
        with open(input["references"], "r") as fh:
            references = set([line.rstrip() for line in fh])

        with open(input["passaged"], "r") as fh:
            with open(output[0], "w") as oh:
                for line in fh:
                    # Parse out strain name from pipe-delimited values.
                    strain, rest = line.split("|", 1)

                    # If this exact strain is not a reference, write it out.
                    if strain not in references:
                        oh.write(line)

rule find_passaged_strains:
    input: "dist/fauna/data/h3n2_ha.fasta"
    output: "outliers/passaged_strains.txt"
    shell: """grep "^>" {input} | grep -v -E "\|(cell|unpassaged)\|" | sed 's/>//' | sort | uniq > {output}"""

rule collect_outliers:
    input:
        "dist/augur/builds/flu/metadata/h3n2_outliers.txt",
        "outliers/non_h3_strains.txt"
    output: "outliers/filtered_strains.txt"
    shell: "sort {input} | uniq > {output}"

rule find_outliers:
    input: "outliers/h3n2_ha.out"
    output: "outliers/non_h3_strains.txt"
    conda: "envs/anaconda.python2.yaml"
    shell: "python scripts/annotate_outliers.py {input} {output}"

rule align_fauna_sequences_to_ncbi_database:
    input:
        fauna_sequences="dist/fauna/data/h3n2_ha.fasta",
        ncbi_database="outliers/iav_ha.nsq"
    output: "outliers/h3n2_ha.out"
    conda: "envs/ncbi.yaml"
    benchmark: "benchmarks/find_outliers.txt"
    threads: 4
    shell: """blastn -num_threads {threads} -task megablast -db outliers/iav_ha -query {input.fauna_sequences} -out {output} -max_target_seqs 5 -outfmt "6 qseqid sseqid pident evalue bitscore stitle" """

rule create_iav_blast_database:
    input: "outliers/iav_ha.fasta"
    output: "outliers/iav_ha.nsq"
    conda: "envs/ncbi.yaml"
    benchmark: "benchmarks/create_iav_blast_database.txt"
    shell: "makeblastdb -dbtype nucl -out iav_ha -in {input}"

rule download_iav_sequences:
    output: "outliers/iav_ha.fasta"
    conda: "envs/ncbi.yaml"
    benchmark: "benchmarks/download_iav_sequences.txt"
    shell: """esearch -db nuccore -query "Influenza A virus[Organism] hemagglutinin" | efetch -format fasta > {output}"""

rule download_sequences_and_titers:
    output:
        "dist/fauna/data/h3n2_ha.fasta",
        "dist/fauna/data/h3n2_na.fasta",
        "dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/download_sequences_and_titers.txt"
    log: "logs/download_sequences_and_titers.log"
    shell: "cd dist/fauna && python download_all.py --virus flu --flu_lineages h3n2 --segments ha na --sequences --titers"

rule clean:
    run:
        print(config.get("clean"))
        if config.get("clean") is not None:
            items = config.get("clean").split(",")
            for year_range in YEAR_RANGES:
                for virus in VIRUSES:
                    for item in items:
                        if item == "all":
                            shell("rm -f dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/processed/flu_h3n2_ha_{year_range}y_{virus}v*".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
                        elif item == "auspice":
                            shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
                        elif item == "tree":
                            shell("rm -f dist/augur/builds/flu/processed/flu_h3n2_ha_{year_range}y_{virus}v*tree*".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/processed/flu_h3n2_ha_{year_range}y_{virus}v*newick*".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
