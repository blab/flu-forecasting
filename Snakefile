# Imports.
import json
import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: clean, download_sequences_and_titers, augur_prepare, summarize_model, aggregate_model_parameters, aggregate_model_accuracy, aggregate_tree_plots

wildcard_constraints:
    sample="\d+",
    viruses="\d+"

# Load configuration parameters.
configfile: "config.json"

SEGMENTS = config["segments"]
YEAR_RANGES = config["year_ranges"]
VIRUSES = config["viruses"]
PREDICTORS = config["predictors"]
NUMBER_OF_SAMPLES = config["number_of_samples"]
SAMPLES = range(NUMBER_OF_SAMPLES)

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
    input: "model_accuracy.tab", "model_parameters.tab", "trees.pdf", "models.tab", "model_fold_change.pdf"

rule trees:
    input: "trees.pdf"

rule aggregate_tree_plots:
    input: expand("figures/trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.pdf", segment=SEGMENTS, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "trees.pdf"
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

rule aggregate_model_fold_change:
    input: expand("model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_fold_change.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule plot_model_fold_change:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    output: "model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python plot_model_fold_change.py {input} {output} {wildcards.year_range} {wildcards.viruses} {wildcards.predictors} {wildcards.sample}"

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
            "clade_error": [accuracy["clade_error"]]
        })
        df.to_csv(output["accuracy"], sep="\t", index=False)

        df = pd.DataFrame(model["params"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample
        df.to_csv(output["parameters"], sep="\t", index=False)

rule run_fitness_model:
    input:
        ha_tree="trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        na_tree="trees/flu_h3n2_na_{year_range}y_{viruses}v_{sample}_tree.json",
        frequencies="frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        titers="dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    output: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    params: predictor_list=_get_predictor_list
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.txt"
    log: "logs/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.log"
    shell: "python fit_model.py {input.ha_tree} {input.frequencies} {output} {params.predictor_list} --na-tree {input.na_tree} --titers {input.titers} &> {log}"

rule estimate_frequencies:
    input:
        tree="trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        weights="region_weights.json"
    output: "frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        start_date=_get_start_date_from_range,
        end_date=_get_end_date_from_range
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/estimate_frequencies_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/estimate_frequencies_{year_range}y_{viruses}v_{sample}.log"
    shell: """python scripts/frequencies.py {input.tree} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--start-date {params.start_date} \
--end-date {params.end_date} \
--weights {input.weights} \
--weights-attribute region \
--include-internal-nodes &> {log}"""

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
        translations = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/translations.json",
        distances = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/distances.json",
        lbi="builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/lbi.json",
        colors = "dist/augur/builds/flu/colors.tsv",
        auspice_config = "config/auspice_config.json"
    output:
        auspice_tree = "trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.json",
        auspice_metadata = "metadata/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_meta.json"
    params:
        geography_traits = "region",
        panels = "tree entropy"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.translations} {input.distances} \
                        {input.lbi} \
            --colors {input.colors} \
            --geography-traits {params.geography_traits} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_metadata} \
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
        tau = 0.2,
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

rule distances:
    message: "Calculating amino acid distances"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        translations = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/translations.json",
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
            {input.translations} \
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

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/tree.nwk",
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/nt_muts.json",
        reference = "dist/augur/builds/flu/metadata/h3n2_{segment}_outgroup.gb"
    output:
        node_data = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aa_muts.json",
        translations = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/translations.json"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
            --alignment-output {output.translations}
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
    threads: 2
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
        sequences = "builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.fasta",
        reference = "dist/augur/builds/flu/metadata/h3n2_{segment}_outgroup.gb"
    output:
        alignment = "builds/results/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}/aligned.fasta"
    conda: "envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps
        """

rule convert_prepared_json_to_metadata_and_sequences:
    input: "dist/augur/builds/flu/prepared/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.json"
    output:
        sequences="builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.fasta",
        metadata="builds/data/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}.tsv"
    conda: "envs/anaconda.python2.yaml"
    log: "logs/convert_prepared_json_{segment}_{year_range}y_{viruses}v_{sample}.log"
    shell: "cd dist/augur/scripts && python prepared_json_to_fasta.py --metadata {SNAKEMAKE_DIR}/{output.metadata} {SNAKEMAKE_DIR}/{input} > {SNAKEMAKE_DIR}/{output.sequences} 2> {SNAKEMAKE_DIR}/{log}"

rule augur_prepare:
    input:
        ha_sequences="dist/fauna/data/h3n2_ha.fasta",
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
  --file_prefix flu_h3n2_*segment*_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample} --lineage h3n2 --segments ha na --time_interval {params.start_date} {params.end_date} --sampling even -r 12y &> {SNAKEMAKE_DIR}/{log}"""

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
