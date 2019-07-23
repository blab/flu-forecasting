"""Rules for generating simulated HA sequences for validation of forecasting models.
"""
BUILD_PATH_SIMULATIONS = "results/builds/{type}/{sample}/"
BUILD_LOG_STEM_SIMULATIONS = "{type}_{sample}"
BUILD_TIMEPOINT_PATH_SIMULATIONS = BUILD_PATH_SIMULATIONS + "timepoints/{timepoint}/"
BUILD_SEGMENT_LOG_STEM_SIMULATIONS = "{type}_{sample}_{timepoint}"


def float_to_datestring(time):
    """Convert a floating point date to a date string

    >>> float_to_datestring(2010.75)
    '2010-10-01'
    >>> float_to_datestring(2011.25)
    '2011-04-01'
    >>> float_to_datestring(2011.0)
    '2011-01-01'
    >>> float_to_datestring(2011.0 + 11.0 / 12)
    '2011-12-01'

    In some cases, the given float value can be truncated leading to unexpected
    conversion between floating point and integer values. This function should
    account for these errors by rounding months to the nearest integer.

    >>> float_to_datestring(2011.9166666666665)
    '2011-12-01'
    >>> float_to_datestring(2016.9609856262834)
    '2016-12-01'
    """
    year = int(time)

    # After accounting for the current year, extract the remainder and convert
    # it to a month using the inverse of the logic used to create the floating
    # point date. If the float date is sufficiently close to the end of the
    # year, rounding can produce a 13th month.
    month = min(int(np.rint(((time - year) * 12) + 1)), 12)

    # Floating point dates do not encode day information, so we always assume
    # they refer to the start of a given month.
    day = 1

    return "%s-%02d-%02d" % (year, month, day)


rule get_strains_by_timepoint_simulated:
    input:
        metadata = _get_metadata_by_wildcards
    output:
        strains = BUILD_TIMEPOINT_PATH_SIMULATIONS + "strains.txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/partition_strains_by_timepoint.py \
            {input.metadata} \
            {wildcards.timepoint} \
            {output}
        """


rule extract_simulated:
    input:
        sequences = _get_sequences_by_wildcards,
        strains = rules.get_strains_by_timepoint_simulated.output.strains,
    output:
        sequences = BUILD_TIMEPOINT_PATH_SIMULATIONS + "filtered_sequences.fasta"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """


rule align_simulated:
    message:
        """
        Aligning sequences for {wildcards}
          - filling gaps with N
        """
    input:
        sequences = rules.extract_simulated.output.sequences,
    output:
        alignment = BUILD_TIMEPOINT_PATH_SIMULATIONS + "aligned.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    threads: 4
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --output {output.alignment} \
            --fill-gaps \
            --nthreads {threads}
        """


rule tree_simulated:
    message: "Building tree ({wildcards})"
    input:
        alignment = rules.align_simulated.output.alignment
    output:
        tree = BUILD_TIMEPOINT_PATH_SIMULATIONS + "tree_raw.nwk"
    conda: "../envs/anaconda.python3.yaml"
    shadow: "minimal"
    benchmark: "benchmarks/tree_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/tree_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method iqtree \
            --nthreads {threads} &> {log}
        """


rule refine_simulated:
    message:
        """
        Refining tree ({wildcards})
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree_simulated.output.tree,
        alignment = rules.align_simulated.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        tree = BUILD_TIMEPOINT_PATH_SIMULATIONS + "tree.nwk",
        node_data = BUILD_TIMEPOINT_PATH_SIMULATIONS + "branch_lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/refine_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/refine_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} &> {log}
        """


rule estimate_frequencies_simulated:
    message:
        """
        Estimating frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree = rules.refine_simulated.output.tree,
        metadata = _get_metadata_by_wildcards,
        weights = "data/region_weights.json"
    output:
        frequencies = BUILD_TIMEPOINT_PATH_SIMULATIONS + "frequencies.json"
    params:
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth = config["frequencies"]["wide_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"],
        pivot_frequency = PIVOT_INTERVAL,
        start_date = _get_start_date_by_wildcards,
        end_date = _get_end_date_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_frequencies_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/estimate_frequencies_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
    shell: """python3 scripts/frequencies.py {input.tree} {input.metadata} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--pivot-frequency {params.pivot_frequency} \
--start-date {params.start_date} \
--end-date {wildcards.timepoint} \
--include-internal-nodes &> {log}"""


rule estimate_diffusion_frequencies_simulated:
    message:
        """
        Estimating diffusion frequencies for {input.tree}
        """
    input:
        tree = rules.refine_simulated.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        frequencies = BUILD_TIMEPOINT_PATH_SIMULATIONS + "diffusion_frequencies.json"
    params:
        pivot_frequency = PIVOT_INTERVAL,
        stiffness = config["frequencies"]["stiffness"],
        inertia = config["frequencies"]["inertia"],
        min_freq = config["frequencies"]["min_freq"],
        min_date = _get_min_date_for_augur_frequencies_by_wildcards,
        max_date = _get_max_date_for_augur_frequencies_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_diffusion_frequencies_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/estimate_diffusion_frequencies_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
    shell: """augur frequencies \
        --method diffusion \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --output {output} \
        --include-internal-nodes \
        --stiffness {params.stiffness} \
        --inertia {params.inertia} \
        --minimal-frequency {params.min_freq} \
        --pivot-interval {params.pivot_frequency} \
        --min-date {params.min_date} \
        --max-date {params.max_date} &> {log}"""


rule ancestral_simulated:
    message: "Reconstructing ancestral sequences and mutations for {wildcards}"
    input:
        tree = rules.refine_simulated.output.tree,
        alignment = rules.align_simulated.output.alignment
    output:
        node_data = BUILD_TIMEPOINT_PATH_SIMULATIONS + "nt_muts.json"
    params:
        inference = "joint"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/ancestral_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/ancestral_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference} &> {log}
        """


rule translate_simulated:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine_simulated.output.tree,
        node_data = rules.ancestral_simulated.output.node_data,
        reference = "config/reference_h3n2_ha.gb"
    output:
        node_data = BUILD_TIMEPOINT_PATH_SIMULATIONS + "aa_muts.json"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/translate_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/translate_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} &> {log}
        """


rule reconstruct_translations_simulated:
    message: "Reconstructing translations for {wildcards.gene}"
    input:
        tree = rules.refine_simulated.output.tree,
        node_data = rules.translate_simulated.output.node_data
    output:
        aa_alignment = BUILD_TIMEPOINT_PATH_SIMULATIONS + "aa-seq_{gene}.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/reconstruct_translations_{gene}_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/reconstruct_translations_{gene}_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes &> {log}
        """


rule convert_translations_to_json_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        translations = translations(segment="ha", path=BUILD_TIMEPOINT_PATH_SIMULATIONS)
    output:
        translations = BUILD_TIMEPOINT_PATH_SIMULATIONS + "aa_seq.json"
    params:
        gene_names = gene_names(segment="ha")
    shell:
        """
        python3 scripts/convert_translations_to_json.py \
            --tree {input.tree} \
            --alignment {input.translations} \
            --gene-names {params.gene_names} \
            --output {output.translations}
        """


rule clades_by_haplotype_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        frequencies = rules.estimate_frequencies_simulated.output.frequencies,
        reference = "config/reference_h3n2_ha.gb",
        translations = translations(segment="ha", path=BUILD_TIMEPOINT_PATH_SIMULATIONS)
    output:
        clades = BUILD_TIMEPOINT_PATH_SIMULATIONS + "clades.json",
        tip_clade_table = BUILD_TIMEPOINT_PATH_SIMULATIONS + "unannotated_tips_to_clades.tsv"
    params:
        gene_names = gene_names(segment="ha"),
        minimum_tips = config["min_tips_per_clade"],
        min_frequency = config["min_frequency_per_clade"]
    conda: "../envs/anaconda.python3.yaml"
    log: "logs/find_clades_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
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


rule delta_frequency_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        frequencies = rules.estimate_diffusion_frequencies_simulated.output.frequencies,
        clades = rules.clades_by_haplotype_simulated.output.clades
    output:
        delta_frequency = BUILD_TIMEPOINT_PATH_SIMULATIONS + "delta_frequency.json"
    params:
        delta_pivots = config["delta_pivots"],
        method = "diffusion"
    conda: "../envs/anaconda.python3.yaml"
    log: "logs/delta_frequency_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
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


rule distances_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        alignments = translations(segment="ha", path=BUILD_TIMEPOINT_PATH_SIMULATIONS),
        distance_maps = _get_distance_maps_for_simulations,
        date_annotations = rules.refine_simulated.output.node_data
    params:
        genes = gene_names(segment="ha"),
        comparisons = _get_distance_comparisons_for_simulations,
        attribute_names = _get_distance_attributes_for_simulations,
        earliest_date = _get_distance_earliest_date_by_wildcards,
        latest_date = _get_distance_latest_date_by_wildcards
    output:
        distances = BUILD_TIMEPOINT_PATH_SIMULATIONS + "distances.json",
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

rule pairwise_distances_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        frequencies = rules.estimate_frequencies_simulated.output.frequencies,
        alignments = translations(segment="ha", path=BUILD_TIMEPOINT_PATH_SIMULATIONS),
        distance_maps = _get_pairwise_distance_maps_for_simulations,
        date_annotations = rules.refine_simulated.output.node_data
    params:
        genes = gene_names(segment="ha"),
        attribute_names = _get_pairwise_distance_attributes_for_simulations,
        years_back_to_compare = config["max_years_for_distances"]
    output:
        distances = BUILD_TIMEPOINT_PATH_SIMULATIONS + "pairwise_distances.json",
    benchmark: "benchmarks/pairwise_distances_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/pairwise_distances_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/pairwise_distances.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --years-back-to-compare {params.years_back_to_compare} \
            --output {output} &> {log}
        """

def _get_cross_immunity_distance_attributes_for_simulations(wildcards):
    return config["cross_immunity"]["h3n2"]["ha"]["distance_attributes"]


def _get_cross_immunity_attributes_for_simulations(wildcards):
    return config["cross_immunity"]["h3n2"]["ha"]["immunity_attributes"]


def _get_cross_immunity_decay_factors_for_simulations(wildcards):
    return config["cross_immunity"]["h3n2"]["ha"]["decay_factors"]


rule cross_immunities_simulated:
    input:
        frequencies = rules.estimate_frequencies_simulated.output.frequencies,
        distances = rules.pairwise_distances_simulated.output.distances
    params:
        distance_attributes = _get_cross_immunity_distance_attributes_for_simulations,
        immunity_attributes = _get_cross_immunity_attributes_for_simulations,
        decay_factors = _get_cross_immunity_decay_factors_for_simulations
    output:
        cross_immunities = BUILD_TIMEPOINT_PATH_SIMULATIONS + "cross_immunity.json",
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


rule lbi_simulated:
    message: "Calculating LBI"
    input:
        tree = rules.refine_simulated.output.tree,
        branch_lengths = rules.refine_simulated.output.node_data
    params:
        tau = config["lbi"]["tau"],
        window = config["lbi"]["window"],
        names = "lbi"
    output:
        lbi = BUILD_TIMEPOINT_PATH_SIMULATIONS + "lbi.json"
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


rule tip_frequencies_simulated:
    message:
        """
        Estimating tip frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree = rules.refine_simulated.output.tree,
        metadata = _get_metadata_by_wildcards,
        weights = "data/region_weights.json"
    output:
        frequencies = "results/auspice/simulated_flu_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + "_tip-frequencies.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=PIVOT_INTERVAL,
        min_date=_get_min_date_for_augur_frequencies_by_wildcards,
        max_date=_get_max_date_for_augur_frequencies_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/tip_frequencies_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/tip_frequencies_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
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


def _get_node_data_for_export_simulated(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by specific builds.
    inputs = [
        rules.refine_simulated.output.node_data,
        rules.ancestral_simulated.output.node_data,
        rules.translate_simulated.output.node_data,
        rules.convert_translations_to_json_simulated.output.translations,
        rules.clades_by_haplotype_simulated.output.clades,
        rules.delta_frequency_simulated.output.delta_frequency,
        rules.distances_simulated.output.distances,
        rules.cross_immunities_simulated.output.cross_immunities,
        rules.lbi_simulated.output.lbi
    ]

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs


rule export_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        metadata = _get_metadata_by_wildcards,
        auspice_config = "config/auspice_config.json",
        node_data = _get_node_data_for_export_simulated,
        colors = "config/colors.tsv"
    output:
        auspice_tree = "results/auspice/simulated_flu_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + "_tree.json",
        auspice_metadata = "results/auspice/simulated_flu_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + "_meta.json"
#        auspice_sequence = "results/auspice/simulated_flu_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + "_seq.json",
    params:
        panels = "tree entropy frequencies"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_metadata} \
            --panels {params.panels} \
            --minify-json
        """


rule convert_node_data_to_table_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        metadata = _get_metadata_by_wildcards,
        node_data = _get_node_data_for_export_simulated
    output:
        table = BUILD_TIMEPOINT_PATH_SIMULATIONS + "node_data.tsv"
    params:
        excluded_fields_arg = _get_excluded_fields_arg
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/node_data_to_table.py \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --jsons {input.node_data} \
            --output {output} \
            {params.excluded_fields_arg} \
            --annotations timepoint={wildcards.timepoint} \
                          lineage=simulated
        """


rule convert_frequencies_to_table_simulated:
    input:
        tree = rules.refine_simulated.output.tree,
        frequencies = rules.estimate_frequencies_simulated.output.frequencies
    output:
        table = BUILD_TIMEPOINT_PATH_SIMULATIONS + "frequencies.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --output {output} \
            --annotations timepoint={wildcards.timepoint}
        """


rule merge_node_data_and_frequencies_simulated:
    input:
        node_data = rules.convert_node_data_to_table_simulated.output.table,
        frequencies = rules.convert_frequencies_to_table_simulated.output.table
    output:
        table = BUILD_TIMEPOINT_PATH_SIMULATIONS + "tip_attributes.tsv"
    run:
        node_data = pd.read_table(input.node_data)
        frequencies = pd.read_table(input.frequencies)
        df = node_data.merge(
            frequencies,
            how="inner",
            on=["strain", "timepoint", "is_terminal"]
        )

        df.to_csv(output.table, sep="\t", index=False, header=True)


def _get_simulated_tip_attributes_by_wildcards(wildcards):
    start_date = config["builds"][wildcards.type][wildcards.sample]["start_date"]
    end_date = config["builds"][wildcards.type][wildcards.sample]["end_date"]

    timepoints_simulations = _get_timepoints_for_build_interval(
        start_date,
        end_date,
        PIVOT_INTERVAL,
        MIN_YEARS_PER_BUILD
    )

    return expand(
        BUILD_PATH_SIMULATIONS.replace("{", "{{").replace("}", "}}") + "timepoints/{timepoint}/tip_attributes.tsv",
        timepoint=timepoints_simulations
    )

rule collect_tip_attributes_simulated:
    input:
        _get_simulated_tip_attributes_by_wildcards
    output:
        attributes = BUILD_PATH_SIMULATIONS + "tip_attributes.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/collect_tables.py \
            --tables {input} \
            --output {output.attributes}
        """


rule target_distances_simulated:
    input:
        attributes = rules.collect_tip_attributes_simulated.output.attributes
    output:
        distances = BUILD_PATH_SIMULATIONS + "target_distances.tsv",
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


rule annotate_naive_tip_attribute_simulated:
    input:
        attributes = rules.collect_tip_attributes_simulated.output.attributes
    output:
        attributes = BUILD_PATH_SIMULATIONS + "tip_attributes_with_naive_predictor.tsv",
    run:
        # Annotate a predictor for a naive model with no growth.
        df = pd.read_csv(input.attributes, sep="\t")
        df["naive"] = 0.0
        df.to_csv(output.attributes, sep="\t", index=False)


rule annotate_weighted_distances_for_tip_attributes_simulated:
    input:
        attributes = rules.annotate_naive_tip_attribute_simulated.output.attributes,
        distances = rules.target_distances_simulated.output.distances
    output:
        attributes = BUILD_PATH_SIMULATIONS + "tip_attributes_with_weighted_distances.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months"]
    shell:
        """
        python3 src/weighted_distances.py \
            --tip-attributes {input.attributes} \
            --distances {input.distances} \
            --delta-months {params.delta_months} \
            --output {output}
        """


# rule collect_annotated_tip_clade_tables_simulated:
#     input:
#         expand("results/builds/{{lineage}}/{{viruses}}_viruses_per_month/{{sample}}/{{start}}--{{end}}/timepoints/{timepoint}/segments/{segment}/tips_to_clades.tsv", timepoint=TIMEPOINTS, segment=SEGMENTS)
#     output:
#         tip_clade_table = BUILD_PATH_SIMULATIONS + "tips_to_clades.tsv"
#     conda: "../envs/anaconda.python3.yaml"
#     shell:
#         """
#         python3 scripts/collect_tables.py \
#             --tables {input} \
#             --output {output.tip_clade_table}
#         """


# rule select_clades_simulated:
#     input:
#         attributes = rules.annotate_naive_tip_attribute_simulated.output.attributes,
#         tips_to_clades = rules.collect_annotated_tip_clade_tables.output.tip_clade_table
#     output:
#         clades = BUILD_PATH_SIMULATIONS + "final_clade_frequencies.tsv"
#     params:
#         primary_segment = config["fitness_model"]["primary_segment"],
#         delta_months = config["fitness_model"]["delta_months"]
#     conda: "../envs/anaconda.python3.yaml"
#     log: "logs/select_clades_" + BUILD_LOG_STEM_SIMULATIONS + ".txt"
#     shell:
#         """
#         python3 scripts/select_clades.py \
#             --tip-attributes {input.attributes} \
#             --tips-to-clades {input.tips_to_clades} \
#             --primary-segment {params.primary_segment} \
#             --delta-months {params.delta_months} \
#             --output {output} &> {log}
#         """


# rule fit_models_by_clades_simulated:
#     input:
#         attributes = rules.annotate_naive_tip_attribute_simulated.output.attributes,
#         final_clade_frequencies = rules.select_clades_simulated.output.clades
#     output:
#         model = BUILD_PATH_SIMULATIONS + "models_by_clades/{predictors}.json"
#     params:
#         predictors = _get_predictor_list,
#         delta_months = config["fitness_model"]["delta_months"],
#         training_window = config["fitness_model"]["training_window"],
#         cost_function = config["fitness_model"]["clade_cost_function"],
#         l1_lambda = config["fitness_model"]["l1_lambda"],
#         pseudocount = config["fitness_model"]["pseudocount"]
#     conda: "../envs/anaconda.python3.yaml"
#     benchmark: "benchmarks/fitness_model_clades_" + BUILD_LOG_STEM_SIMULATIONS + "_{predictors}.txt"
#     log: "logs/fitness_model_clades_" + BUILD_LOG_STEM_SIMULATIONS + "_{predictors}.txt"
#     shell:
#         """
#         python3 src/fit_model.py \
#             --tip-attributes {input.attributes} \
#             --final-clade-frequencies {input.final_clade_frequencies} \
#             --training-window {params.training_window} \
#             --delta-months {params.delta_months} \
#             --predictors {params.predictors} \
#             --cost-function {params.cost_function} \
#             --l1-lambda {params.l1_lambda} \
#             --pseudocount {params.pseudocount} \
#             --target clades \
#             --output {output} &> {log}
#         """


rule fit_models_by_distances_simulated:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes_simulated.output.attributes,
        distances = rules.target_distances_simulated.output.distances
    output:
        model = BUILD_PATH_SIMULATIONS + "models_by_distances/{predictors}.json"
    params:
        predictors = _get_predictor_list,
        delta_months = config["fitness_model"]["delta_months"],
        training_window = config["fitness_model"]["training_window"],
        cost_function = config["fitness_model"]["distance_cost_function"],
        l1_lambda = config["fitness_model"]["l1_lambda"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/fitness_model_distances_" + BUILD_LOG_STEM_SIMULATIONS + "_{predictors}.txt"
    log: "logs/fitness_model_distances_" + BUILD_LOG_STEM_SIMULATIONS + "_{predictors}.txt"
    shell:
        """
        python3 src/fit_model.py \
            --tip-attributes {input.attributes} \
            --training-window {params.training_window} \
            --delta-months {params.delta_months} \
            --predictors {params.predictors} \
            --cost-function {params.cost_function} \
            --l1-lambda {params.l1_lambda} \
            --target distances \
            --distances {input.distances} \
            --output {output} &> {log}
        """


rule plot_tree_simulated:
    input:
        auspice_tree = rules.export_simulated.output.auspice_tree
    output:
        tree = "results/figures/trees/flu_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + "_tree.pdf"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/plot_tree_simulated_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".txt"
    log: "logs/plot_tree_simulated_" + BUILD_SEGMENT_LOG_STEM_SIMULATIONS + ".log"
    params:
        start = _get_start_date_by_wildcards,
        end = _get_end_date_by_wildcards
    shell:
        """
        python3 scripts/plot_tree.py \
            {input} \
            {output} \
            --start-date {params.start} \
            --end-date {params.end} &> {log}
        """

# rule aggregate_tree_plots_simulated:
#     input: expand(rules.plot_tree_simulated.output.tree, percentage=PERCENTAGE, start=START_DATE_SIMULATIONS, end=END_DATE_SIMULATIONS, timepoint=TIMEPOINTS_SIMULATIONS)
#     output:
#         trees="results/figures/trees_simulated.pdf"
#     shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"
