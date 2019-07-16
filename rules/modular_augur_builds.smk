"""
Rules to build auspice JSONs from sequences and titers using modular augur.
"""
# Imports.
import json
import pandas as pd
from pathlib import Path
from snakemake.exceptions import IncompleteCheckpointException


rule download_sequences:
    message: "Downloading {wildcards.lineage} {wildcards.segment} sequences from fauna"
    output:
        sequences = "data/{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/download_sequences_{lineage}_{segment}.txt"
    log: "logs/download_sequences_{lineage}_{segment}.log"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            --path data \
            --fstem {wildcards.lineage}_{wildcards.segment}
        """

rule download_all_titers_by_assay:
    message: "Downloading {wildcards.lineage} {wildcards.assay} titers from fauna"
    output:
        titers = "data/{lineage}_{assay}_titers.json"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/download_all_titers_{lineage}_{assay}.txt"
    log: "logs/download_all_titers_{lineage}_{assay}.log"
    params:
        databases = "tdb cdc_tdb crick_tdb vidrl_tdb niid_tdb"
    shell:
        """
        python3 {path_to_fauna}/tdb/download.py \
            --database {params.databases} \
            --virus flu \
            --subtype {wildcards.lineage} \
            --select assay_type:{wildcards.assay} \
            --path data \
            --fstem {wildcards.lineage}_{wildcards.assay} \
            --ftype json
        """

rule get_titers_by_passage:
    message: "Parsing {wildcards.passage}-passaged titers for {wildcards.lineage} {wildcards.assay}"
    input:
        titers = rules.download_all_titers_by_assay.output.titers
    output:
        titers = "data/{lineage}_{passage}_{assay}_titers.tsv"
    benchmark: "benchmarks/get_titers_{lineage}_{passage}_{assay}.txt"
    log: "logs/get_titers_{lineage}_{passage}_{assay}.log"
    run:
        df = pd.read_json(input.titers)
        passaged = (df["serum_passage_category"] == wildcards.passage)
        tdb_passaged = df["index"].apply(lambda index: isinstance(index, list) and wildcards.passage in index)
        tsv_fields = [
            "virus_strain",
            "serum_strain",
            "serum_id",
            "source",
            "titer",
            "assay_type"
        ]

        titers_df = df.loc[(passaged | tdb_passaged), tsv_fields]
        titers_df.to_csv(output.titers, sep="\t", header=False, index=False)

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download_sequences.output.sequences
    output:
        sequences = "results/builds/sequences_{lineage}_{segment}.fasta",
        metadata = "results/builds/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division location passage authors age gender"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message:
        """
        Filtering {wildcards.lineage} {wildcards.segment} sequences:
          - less than {params.min_length} bases
          - outliers
          - samples with missing region and country metadata
        """
    input:
        metadata = rules.parse.output.metadata,
        sequences = rules.parse.output.sequences,
        exclude = "config/outliers_{lineage}.txt"
    output:
        sequences = 'results/builds/filtered_{lineage}_{segment}.fasta'
    params:
        min_length = MIN_LENGTH
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/filter_{lineage}_{segment}.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --exclude {input.exclude} \
            --exclude-where country=? region=? passage=egg \
            --output {output}
        """

rule filter_metadata:
    message:
        """
        Excluding strains with ambiguous dates for {wildcards.lineage} {wildcards.segment}
        """
    input:
        metadata = rules.parse.output.metadata,
    output:
        metadata = 'results/builds/filtered_metadata_{lineage}_{segment}.tsv'
    run:
        df = pd.read_csv(input.metadata, sep="\t")
        df[~df["date"].str.contains("XX")].to_csv(output.metadata, sep="\t", header=True, index=False)

rule select_strains:
    input:
        sequences = expand("results/builds/filtered_{{lineage}}_{segment}.fasta", segment=SEGMENTS),
        metadata = expand("results/builds/filtered_metadata_{{lineage}}_{segment}.tsv", segment=SEGMENTS),
        titers = expand("data/{{lineage}}_{passage}_{assay}_titers.tsv", passage=TITER_PASSAGES, assay=TITER_ASSAYS),
        include = "config/references_{lineage}.txt"
    output:
        strains = "results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/strains.txt"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/select_strains_{lineage}_{viruses}v_{sample}_{start}_{end}.txt"
    log: "logs/select_strains_{lineage}_{viruses}v_{sample}_{start}_{end}.log"
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {SEGMENTS} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --time-interval {wildcards.start} {wildcards.end} \
            --viruses_per_month {wildcards.viruses} \
            --titers {input.titers} \
            --output {output.strains}
        """

rule extract_strain_metadata:
    input:
        strains = rules.select_strains.output.strains,
        metadata = "results/builds/metadata_{lineage}_ha.tsv"
    output:
        metadata = "results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/strains_metadata.tsv"
    run:
        strains = pd.read_table(input.strains, header=None, names=["strain"])
        metadata = pd.read_table(input.metadata)
        selected_metadata = strains.merge(metadata, how="left", on="strain")
        selected_metadata.to_csv(output.metadata, sep="\t", index=False)

rule get_strains_by_timepoint:
    input:
        metadata = rules.extract_strain_metadata.output.metadata
    output:
        strains = "results/builds/{lineage}/{viruses}_viruses_per_month/{sample}/{start}--{end}/timepoints/{timepoint}/strains.txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/partition_strains_by_timepoint.py \
            {input.metadata} \
            {wildcards.timepoint} \
            {output}
        """

rule extract:
    input:
        sequences = rules.filter.output.sequences,
        strains = rules.get_strains_by_timepoint.output.strains
    output:
        sequences = BUILD_SEGMENT_PATH + "filtered_sequences.fasta"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference} for {wildcards}
          - filling gaps with N
        """
    input:
        sequences = rules.extract.output.sequences,
        reference = "config/reference_{lineage}_{segment}.gb"
    output:
        alignment = BUILD_SEGMENT_PATH + "aligned.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_" + BUILD_SEGMENT_LOG_STEM + ".txt"
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

rule tree:
    message: "Building tree ({wildcards})"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = BUILD_SEGMENT_PATH + "tree_raw.nwk"
    conda: "../envs/anaconda.python3.yaml"
    shadow: "minimal"
    benchmark: "benchmarks/tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --nthreads {threads} &> {log}
        """

rule refine:
    message:
        """
        Refining tree ({wildcards})
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - use fixed clock rate of {params.clock_rate}
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.parse.output.metadata
    output:
        tree = BUILD_SEGMENT_PATH + "tree.nwk",
        node_data = BUILD_SEGMENT_PATH + "branch_lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = _get_clock_rate_by_wildcards,
        clock_std_dev = _get_clock_std_dev_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/refine_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/refine_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} &> {log}
        """

rule estimate_frequencies:
    message:
        """
        Estimating frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree=rules.refine.output.tree,
        metadata=rules.parse.output.metadata,
        weights="data/region_weights.json"
    output:
        frequencies = BUILD_SEGMENT_PATH + "frequencies.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=PIVOT_INTERVAL
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/estimate_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell: """python3 scripts/frequencies.py {input.tree} {input.metadata} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--pivot-frequency {params.pivot_frequency} \
--start-date {wildcards.start} \
--end-date {wildcards.timepoint} \
--include-internal-nodes &> {log}"""
#--weights {input.weights} \
#--weights-attribute region \

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations for {wildcards}"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = BUILD_SEGMENT_PATH + "nt_muts.json"
    params:
        inference = "joint"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/ancestral_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/ancestral_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference} &> {log}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "config/reference_{lineage}_{segment}.gb"
    output:
        node_data = BUILD_SEGMENT_PATH + "aa_muts.json"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/translate_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/translate_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} &> {log}
        """

rule reconstruct_translations:
    message: "Reconstructing translations for {wildcards.gene}"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.translate.output.node_data
    output:
        aa_alignment = BUILD_SEGMENT_PATH + "aa-seq_{gene}.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/reconstruct_translations_{gene}_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/reconstruct_translations_{gene}_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes &> {log}
        """

rule convert_translations_to_json:
    input:
        tree = rules.refine.output.tree,
        translations = translations
    output:
        translations = BUILD_SEGMENT_PATH + "aa_seq.json"
    params:
        gene_names = gene_names
    shell:
        """
        python3 scripts/convert_translations_to_json.py \
            --tree {input.tree} \
            --alignment {input.translations} \
            --gene-names {params.gene_names} \
            --output {output.translations}
        """

rule mutation_frequencies:
    message:
        """
        Estimating diffusion-based mutation frequencies for {wildcards.lineage} {wildcards.segment} genes: {params.gene_names}
        """
    input:
        metadata = rules.parse.output.metadata,
        translations = translations
    output:
        frequencies = BUILD_SEGMENT_PATH + "mutation_frequencies.json"
    params:
        gene_names = gene_names,
        pivot_frequency = PIVOT_INTERVAL,
        min_date = _get_min_date_for_augur_frequencies,
        max_date = _get_max_date_for_augur_frequencies,
        min_frequency = config["frequencies"]["min_mutation_frequency"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/mutation_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/mutation_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        augur frequencies \
            --method diffusion \
            --metadata {input.metadata} \
            --alignments {input.translations} \
            --gene-names {params.gene_names} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --minimal-frequency {params.min_frequency} \
            --pivot-interval {params.pivot_frequency} \
            --output {output}
        """

rule clades_by_haplotype:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_frequencies.output.frequencies,
        reference = "config/reference_{lineage}_{segment}.gb",
        translations = translations
    output:
        clades = BUILD_SEGMENT_PATH + "clades.json",
        tip_clade_table = BUILD_SEGMENT_PATH + "unannotated_tips_to_clades.tsv"
    params:
        gene_names = gene_names,
        minimum_tips = config["min_tips_per_clade"],
        min_frequency = config["min_frequency_per_clade"]
    conda: "../envs/anaconda.python3.yaml"
    log: "logs/find_clades_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        python3 scripts/find_clades.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --reference {input.reference} \
            --translations {input.translations} \
            --gene-names {params.gene_names} \
            --minimum-tips {params.minimum_tips} \
            --minimum-frequency {params.min_frequency} \
            --use-hash-ids \
            --output {output.clades} \
            --output-tip-clade-table {output.tip_clade_table} &> {log}
        """

rule annotate_tip_clade_table:
    input:
        tip_clade_table = rules.clades_by_haplotype.output.tip_clade_table
    output:
        tip_clade_table = BUILD_SEGMENT_PATH + "tips_to_clades.tsv"
    run:
        df = pd.read_table(input["tip_clade_table"])
        df["lineage"] = wildcards.lineage
        df["segment"] = wildcards.segment
        df["timepoint"] = wildcards.timepoint
        df.to_csv(output["tip_clade_table"], sep="\t", index=False)

rule estimate_diffusion_frequencies:
    message:
        """
        Estimating diffusion frequencies for {input.tree}
        """
    input:
        tree=rules.refine.output.tree,
        metadata=rules.parse.output.metadata
    output:
        frequencies = BUILD_SEGMENT_PATH + "diffusion_frequencies.json"
    params:
        pivot_frequency = PIVOT_INTERVAL,
        stiffness = config["frequencies"]["stiffness"],
        inertia = config["frequencies"]["inertia"],
        min_freq = config["frequencies"]["min_freq"],
        min_date = _get_min_date_for_augur_frequencies,
        max_date = _get_max_date_for_augur_frequencies
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_diffusion_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/estimate_diffusion_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell: """augur frequencies \
        --method diffusion \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --output {output} \
        --include-internal-nodes \
        --ignore-char X \
        --stiffness {params.stiffness} \
        --inertia {params.inertia} \
        --minimal-frequency {params.min_freq} \
        --pivot-interval {params.pivot_frequency} \
        --min-date {params.min_date} \
        --max-date {params.max_date} &> {log}"""

rule delta_frequency:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_diffusion_frequencies.output.frequencies,
        clades = rules.clades_by_haplotype.output.clades
    output:
        delta_frequency = BUILD_SEGMENT_PATH + "delta_frequency.json"
    params:
        delta_pivots = config["delta_pivots"],
        method = "diffusion"
    conda: "../envs/anaconda.python3.yaml"
    log: "logs/delta_frequency_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        python3 scripts/calculate_delta_frequency.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --frequency-method {params.method} \
            --clades {input.clades} \
            --delta-pivots {params.delta_pivots} \
            --output {output.delta_frequency} &> {log}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = BUILD_SEGMENT_PATH + "traits.json",
    params:
        columns = "region country"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/traits_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = _get_distance_maps_by_lineage_and_segment,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names,
        comparisons = _get_distance_comparisons_by_lineage_and_segment,
        attribute_names = _get_distance_attributes_by_lineage_and_segment,
        earliest_date = _get_distance_earliest_date_by_wildcards,
        latest_date = _get_distance_latest_date_by_wildcards
    output:
        distances = BUILD_SEGMENT_PATH + "distances.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --earliest-date {params.earliest_date} \
            --latest-date {params.latest_date} \
            --output {output}
        """

def _get_cross_immunity_distance_attributes_by_lineage_and_segment(wildcards):
    return config["cross_immunity"][wildcards.lineage][wildcards.segment]["distance_attributes"]

def _get_cross_immunity_attributes_by_lineage_and_segment(wildcards):
    return config["cross_immunity"][wildcards.lineage][wildcards.segment]["immunity_attributes"]

def _get_cross_immunity_decay_factors_by_lineage_and_segment(wildcards):
    return config["cross_immunity"][wildcards.lineage][wildcards.segment]["decay_factors"]

rule cross_immunities:
    input:
        frequencies = rules.estimate_frequencies.output.frequencies,
        distances = rules.distances.output.distances
    params:
        distance_attributes = _get_cross_immunity_distance_attributes_by_lineage_and_segment,
        immunity_attributes = _get_cross_immunity_attributes_by_lineage_and_segment,
        decay_factors = _get_cross_immunity_decay_factors_by_lineage_and_segment
    output:
        cross_immunities = BUILD_SEGMENT_PATH + "cross_immunity.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/cross_immunity.py \
            --frequencies {input.frequencies} \
            --distances {input.distances} \
            --distance-attributes {params.distance_attributes} \
            --immunity-attributes {params.immunity_attributes} \
            --decay-factors {params.decay_factors} \
            --output {output}
        """

rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = config["lbi"]["tau"],
        window = config["lbi"]["window"],
        names = "lbi"
    output:
        lbi = BUILD_SEGMENT_PATH + "lbi.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window}
        """

def _get_min_date_for_translation_filter(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    min_date = timepoint - pd.DateOffset(years=config["years_for_titer_alignments"])
    return min_date.strftime("%Y-%m-%d")

rule filter_translations_by_date:
    input:
        alignments = rules.reconstruct_translations.output.aa_alignment,
        branch_lengths = rules.refine.output.node_data
    output:
        alignments = BUILD_SEGMENT_PATH + "filtered-aa-seq_{gene}.fasta"
    params:
        min_date = _get_min_date_for_translation_filter
    shell:
        """
        python3 scripts/filter_translations.py \
            --alignment {input.alignments} \
            --branch-lengths {input.branch_lengths} \
            --min-date {params.min_date} \
            --output {output}
        """

rule titers_sub:
    input:
        titers = expand("data/{{lineage}}_{passage}_{assay}_titers.tsv", passage=TITER_PASSAGES, assay=TITER_ASSAYS),
        aa_muts = rules.translate.output,
        alignments = filtered_translations,
        tree = rules.refine.output.tree
    params:
        genes = gene_names
    output:
        titers_model = BUILD_SEGMENT_PATH + "titers-sub-model.json",
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/titers_sub_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/titers_sub_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --allow-empty-model \
            --output {output.titers_model} &> {log}
        """

rule titers_tree:
    input:
        titers = expand("data/{{lineage}}_{passage}_{assay}_titers.tsv", passage=TITER_PASSAGES, assay=TITER_ASSAYS),
        tree = rules.refine.output.tree
    output:
        titers_model = BUILD_SEGMENT_PATH + "titers-tree-model.json",
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/titers_tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/titers_tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --allow-empty-model \
            --output {output.titers_model} &> {log}
        """

rule convert_titer_model_to_distance_map:
    input:
        model = rules.titers_sub.output.titers_model
    output:
        distance_map = BUILD_SEGMENT_PATH + "titer_substitution_distance_map.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/titer_model_to_distance_map.py \
            --model {input.model} \
            --output {output}
        """

rule titer_distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        distance_maps = rules.convert_titer_model_to_distance_map.output.distance_map,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names,
        comparisons = "root ancestor pairwise",
        attribute_names = "cTiterSub cTiterSub_star cTiterSub_pairwise",
        earliest_date = _get_distance_earliest_date_by_wildcards,
        latest_date = _get_distance_latest_date_by_wildcards
    output:
        distances = BUILD_SEGMENT_PATH + "titer_substitution_distances.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} {input.distance_maps} {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --earliest-date {params.earliest_date} \
            --latest-date {params.latest_date} \
            --output {output}
        """

rule titer_cross_immunities:
    input:
        frequencies = rules.estimate_frequencies.output.frequencies,
        distances = rules.titer_distances.output.distances
    params:
        distance_attributes = "cTiterSub_pairwise",
        immunity_attributes = "cTiterSub_x",
        decay_factors = "14.0"
    output:
        cross_immunities = BUILD_SEGMENT_PATH + "titer_substitution_cross_immunity.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/cross_immunity.py \
            --frequencies {input.frequencies} \
            --distances {input.distances} \
            --distance-attributes {params.distance_attributes} \
            --immunity-attributes {params.immunity_attributes} \
            --decay-factors {params.decay_factors} \
            --output {output}
        """

rule tip_frequencies:
    message:
        """
        Estimating tip frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree=rules.refine.output.tree,
        metadata=rules.parse.output.metadata,
        weights="data/region_weights.json"
    output:
        frequencies = "results/auspice/flu_{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}_tip-frequencies.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=PIVOT_INTERVAL,
        min_date=_get_min_date_for_augur_frequencies,
        max_date=_get_max_date_for_augur_frequencies
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/tip_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/tip_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivot-interval {params.pivot_frequency} \
            --output {output}
        """
#            --weights {input.weights} \
#            --weights-attribute region \

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by specific builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.convert_translations_to_json.output.translations,
        rules.traits.output.node_data,
        rules.clades_by_haplotype.output.clades,
        rules.lbi.output.lbi,
        rules.delta_frequency.output.delta_frequency
    ]

    # Define segment-specific inputs for auspice JSONs.
    # For example, antigenic assays only make sense for HA and NA.
    if wildcards.lineage in ["h3n2"] and wildcards.segment in ["ha", "na"]:
        inputs.extend([
            rules.titers_tree.output.titers_model,
            rules.titers_sub.output.titers_model,
            rules.cross_immunities.output.cross_immunities,
            rules.titer_distances.output.distances,
            rules.titer_cross_immunities.output.cross_immunities
        ])

    # If the current lineage and segment have a mask configuration defined,
    # calculate distances.
    build_mask_config = _get_build_mask_config(wildcards)
    if build_mask_config is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = "config/auspice_config.json",
        node_data = _get_node_data_for_export,
        colors = "config/colors.tsv"
    output:
        auspice_tree = "results/auspice/flu_{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}_tree.json",
        auspice_metadata = "results/auspice/flu_{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}_meta.json",
        auspice_sequence = "results/auspice/flu_{lineage}_{viruses}_{sample}_{start}_{end}_{timepoint}_{segment}_seq.json",
    params:
        geography_traits = "region",
        panels = "tree entropy"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --geography-traits {params.geography_traits} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_metadata} \
            --output-sequence {output.auspice_sequence} \
            --panels {params.panels} \
            --minify-json
        """

rule convert_node_data_to_table:
    input:
        tree = rules.refine.output.tree,
        node_data = _get_node_data_for_export
    output:
        table = BUILD_SEGMENT_PATH + "node_data.tsv"
    params:
        excluded_fields_arg = _get_excluded_fields_arg
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/node_data_to_table.py \
            --tree {input.tree} \
            --jsons {input.node_data} \
            --output {output} \
            {params.excluded_fields_arg} \
            --annotations timepoint={wildcards.timepoint} \
                          lineage={wildcards.lineage} \
                          segment={wildcards.segment}
        """

rule convert_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_frequencies.output.frequencies
    output:
        table = BUILD_SEGMENT_PATH + "frequencies.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --output {output} \
            --annotations timepoint={wildcards.timepoint}
        """

rule merge_node_data_and_frequencies:
    input:
        node_data = rules.convert_node_data_to_table.output.table,
        frequencies = rules.convert_frequencies_to_table.output.table
    output:
        table = BUILD_SEGMENT_PATH + "tip_attributes.tsv"
    run:
        node_data = pd.read_table(input.node_data)
        frequencies = pd.read_table(input.frequencies)
        df = node_data.merge(
            frequencies,
            how="inner",
            on=["strain", "timepoint", "is_terminal"]
        )

        df.to_csv(output.table, sep="\t", index=False, header=True)

rule collect_tip_attributes:
    input:
        expand("results/builds/{{lineage}}/{{viruses}}_viruses_per_month/{{sample}}/{{start}}--{{end}}/timepoints/{timepoint}/segments/{segment}/tip_attributes.tsv", timepoint=TIMEPOINTS, segment=SEGMENTS)
    output:
        attributes = BUILD_PATH + "tip_attributes.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/collect_tables.py \
            --tables {input} \
            --output {output.attributes}
        """

rule collect_annotated_tip_clade_tables:
    input:
        expand("results/builds/{{lineage}}/{{viruses}}_viruses_per_month/{{sample}}/{{start}}--{{end}}/timepoints/{timepoint}/segments/{segment}/tips_to_clades.tsv", timepoint=TIMEPOINTS, segment=SEGMENTS)
    output:
        tip_clade_table = BUILD_PATH + "tips_to_clades.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/collect_tables.py \
            --tables {input} \
            --output {output.tip_clade_table}
        """

rule target_distances:
    input:
        attributes = rules.collect_tip_attributes.output.attributes
    output:
        distances = BUILD_PATH + "target_distances.tsv",
    params:
        delta_months = config["fitness_model"]["delta_months"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --output {output}
        """
