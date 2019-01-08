"""
Rules to build auspice JSONs from sequences and titers using modular augur.
"""

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
            --database cdc_tdb \
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
            --exclude-where country=? region=? \
            --output {output}
        """

rule select_strains:
    input:
        sequences = expand("results/builds/filtered_{{lineage}}_{segment}.fasta", segment=SEGMENTS),
        metadata = expand("results/builds/metadata_{{lineage}}_{segment}.tsv", segment=SEGMENTS),
        titers = rules.download_titers.output.titers,
        include = "config/references_{lineage}.txt"
    output:
        strains = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/strains.txt",
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/select_strains_{lineage}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/select_strains_{lineage}_{year_range}y_{viruses}v_{sample}.log"
    params:
        start_date=_get_start_date_from_range,
        end_date=_get_end_date_from_range
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {SEGMENTS} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --time-interval {params.start_date} {params.end_date} \
            --viruses_per_month {wildcards.viruses} \
            --titers {input.titers} \
            --output {output.strains}
        """

rule extract:
    input:
        sequences = rules.filter.output.sequences,
        strains = rules.select_strains.output.strains
    output:
        sequences = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/filtered_sequences.fasta"
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
        Aligning sequences to {input.reference} for {wildcards.segment}_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample}
          - filling gaps with N
        """
    input:
        sequences = rules.extract.output.sequences,
        reference = "config/{lineage}_{segment}_outgroup.gb"
    output:
        alignment = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/aligned.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
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
    message: "Building tree ({wildcards.lineage}_{wildcards.segment}_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample})"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/tree_raw.nwk"
    conda: "../envs/anaconda.python3.yaml"
    shadow: "minimal"
    benchmark: "benchmarks/tree_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
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
        Refining tree ({wildcards.lineage}_{wildcards.segment}_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample})
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
        tree = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/tree.nwk",
        node_data = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/branch_lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = _get_clock_rate_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/refine_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
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

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/nt_muts.json"
    params:
        inference = "joint"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/ancestral_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
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
        node_data = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/aa_muts.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/builds/flu_{lineage}_{year_range}y_{viruses}v_{sample}/{segment}/traits.json",
    params:
        columns = "region country"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/traits_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.traits.output.node_data
        # Omit these annotations for now
        # rules.titers_tree.output.titers_model,
        # rules.titers_sub.output.titers_model,
        # rules.clades.output.clades,
        # rules.lbi.output.lbi
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
        auspice_tree = "results/auspice/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}_tree.json",
        auspice_metadata = "results/auspice/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}_meta.json",
        auspice_sequence = "results/auspice/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}_seq.json",
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
