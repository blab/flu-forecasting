"""
Rules to build auspice JSONs from sequences and titers using modular augur.
"""
# Imports.
from augur.frequency_estimators import get_pivots, timestamp_to_float
import json
import pandas as pd
from pathlib import Path
from snakemake.exceptions import IncompleteCheckpointException


rule download_sequences:
    message: "Downloading {wildcards.segment} sequences from fauna"
    output:
        sequences = "data/h3n2_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/download_sequences_{segment}.txt"
    log: "logs/download_sequences_{segment}.log"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_h3n2 \
            --path data \
            --fstem h3n2_{wildcards.segment}
        """

rule download_titers:
    message: "Downloading titers from fauna"
    output:
        titers = "data/h3n2_hi_titers.tsv"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/download_titers.txt"
    log: "logs/download_titers.log"
    shell:
        """
        python3 {path_to_fauna}/tdb/download.py \
            --database tdb cdc_tdb \
            --virus flu \
            --subtype h3n2 \
            --select assay_type:hi \
            --path data \
            --fstem h3n2_hi
        """

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

rule select_strains:
    input:
        sequences = expand("results/builds/filtered_{{lineage}}_{segment}.fasta", segment=SEGMENTS),
        metadata = expand("results/builds/metadata_{{lineage}}_{segment}.tsv", segment=SEGMENTS),
        titers = rules.download_titers.output.titers,
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
        reference = "config/{lineage}_{segment}_outgroup.gb"
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
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --nthreads {threads}
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
        clock_rate = _get_clock_rate_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/refine_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
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
--weights {input.weights} \
--weights-attribute region \
--include-internal-nodes &> {log}"""

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
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "config/{lineage}_{segment}_outgroup.gb"
    output:
        node_data = BUILD_SEGMENT_PATH + "aa_muts.json"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/translate_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule reconstruct_translations:
    message: "Reconstructing translations for {wildcards.gene}"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.translate.output.node_data
    output:
        aa_alignment = BUILD_SEGMENT_PATH + "aa-seq_{gene}.fasta"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes
        """

genes_to_translate = {'ha':['SigPep', 'HA1', 'HA2'], 'na':['NA']}
def gene_names(wildcards):
    return genes_to_translate[wildcards.segment]

def translations(wildcards):
    genes = gene_names(wildcards)
    return [BUILD_SEGMENT_PATH + "aa-seq_%s.fasta" % gene
            for gene in genes]

rule clades_by_haplotype:
    input:
        tree = rules.refine.output.tree,
        translations = translations,
    output:
        clades = BUILD_SEGMENT_PATH + "clades.json"
    params:
        gene_names = gene_names,
        minimum_tips = config["min_tips_per_clade"]
    shell:
        """
        python3 scripts/find_clades.py \
            --tree {input.tree} \
            --translations {input.translations} \
            --gene-names {params.gene_names} \
            --minimum-tips {params.minimum_tips} \
            --output {output.clades}
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
        masks = "config/{segment}_masks.tsv"
    params:
        genes = gene_names,
        attribute_names = _get_mask_attribute_names_by_wildcards,
        mask_names = _get_mask_names_by_wildcards
    output:
        distances = BUILD_SEGMENT_PATH + "distances.json",
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --masks {input.masks} \
            --output {output} \
            --attribute-names {params.attribute_names} \
            --mask-names {params.mask_names}
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

rule titers_sub:
    input:
        titers = rules.download_titers.output.titers,
        aa_muts = rules.translate.output,
        alignments = translations,
        tree = rules.refine.output.tree
    params:
        genes = gene_names
    output:
        titers_model = BUILD_SEGMENT_PATH + "titers-sub-model.json",
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule titers_tree:
    input:
        titers = rules.download_titers.output.titers,
        tree = rules.refine.output.tree
    output:
        titers_model = BUILD_SEGMENT_PATH + "titers-tree-model.json",
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.traits.output.node_data,
        rules.clades_by_haplotype.output.clades,
        rules.lbi.output.lbi,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model
    ]

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
            --panels {params.panels}
        """

rule convert_node_data_to_table:
    input:
        tree = rules.refine.output.tree,
        node_data = [
            rules.lbi.output.lbi,
            rules.distances.output.distances,
            rules.titers_tree.output.titers_model,
            rules.titers_sub.output.titers_model
        ]
    output:
        table = BUILD_SEGMENT_PATH + "node_data.tsv"
    shell:
        """
        python3 scripts/node_data_to_table.py \
            --tree {input.tree} \
            --jsons {input.node_data} \
            --output {output} \
            --annotations timepoint={wildcards.timepoint}
        """

rule convert_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_frequencies.output.frequencies
    output:
        table = BUILD_SEGMENT_PATH + "frequencies.tsv"
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
    run:
        # Concatenate tip attributes across all timepoints.
        df = pd.concat([pd.read_table(i) for i in input], ignore_index=True)

        # Filter out any tips that have a frequency of zero.
        df = df[df["frequency"] > 0].copy()

        df.to_csv(output["attributes"], sep="\t", index=False)
