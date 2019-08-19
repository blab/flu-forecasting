"""Rules for generating simulated HA sequences for validation of forecasting models.
"""
BUILD_PATH = "results/builds/{type}/{sample}/"
BUILD_LOG_STEM = "{type}_{sample}"
BUILD_TIMEPOINT_PATH = BUILD_PATH + "timepoints/{timepoint}/"
BUILD_SEGMENT_LOG_STEM = "{type}_{sample}_{timepoint}"


rule get_titer_sequences_by_timepoint:
    input:
        titers = _get_titers_by_wildcards,
        metadata = _get_complete_metadata_by_wildcards,
        sequences = _get_sequences_by_wildcards
    output:
        sequences = BUILD_TIMEPOINT_PATH + "titer_sequences.fasta"
    params:
        years_back = config["years_for_titer_alignments"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/select_titer_strain_sequences.py \
            --titers {input.titers} \
            --metadata {input.metadata} \
            --sequences {input.sequences} \
            --timepoint {wildcards.timepoint} \
            --years-back {params.years_back} \
            --output {output.sequences}
        """


rule align_titer_sequences:
    input:
        sequences = rules.get_titer_sequences_by_timepoint.output.sequences,
        reference = _get_reference
    output:
        alignment = BUILD_TIMEPOINT_PATH + "aligned_titer_sequences.fasta"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/align_titer_sequences" + BUILD_SEGMENT_LOG_STEM + ".txt"
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


rule translate_titer_sequences:
    input:
        alignment = rules.align_titer_sequences.output.alignment,
        reference = _get_reference
    output:
        sequences = BUILD_TIMEPOINT_PATH + "titer-sequence-aa-seq_{gene}.fasta"
    shell:
        """
        python3 scripts/translate.py \
            --sequences {input.alignment} \
            --reference-sequence {input.reference} \
            --genes {wildcards.gene} \
            --output {output.sequences}
        """


rule get_strains_by_timepoint:
    input:
        metadata = _get_metadata_by_wildcards
    output:
        strains = BUILD_TIMEPOINT_PATH + "strains.txt"
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
        sequences = _get_sequences_by_wildcards,
        strains = rules.get_strains_by_timepoint.output.strains,
    output:
        sequences = BUILD_TIMEPOINT_PATH + "filtered_sequences.fasta"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """


rule align:
    input:
        sequences = rules.extract.output.sequences,
        reference = _get_reference
    output:
        alignment = BUILD_TIMEPOINT_PATH + "aligned.fasta"
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
        tree = BUILD_TIMEPOINT_PATH + "tree_raw.nwk"
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
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = _get_metadata_by_wildcards
    output:
        tree = BUILD_TIMEPOINT_PATH + "tree.nwk",
        node_data = BUILD_TIMEPOINT_PATH + "branch_lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_rate = _get_clock_rate_argument,
        clock_std_dev = _get_clock_std_dev_argument
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
            {params.clock_rate} \
            {params.clock_std_dev} \
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
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        weights = "data/region_weights.json"
    output:
        frequencies = BUILD_TIMEPOINT_PATH + "frequencies.json"
    params:
        narrow_bandwidth = config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth = config["frequencies"]["wide_bandwidth"],
        proportion_wide = config["frequencies"]["proportion_wide"],
        pivot_frequency = _get_pivot_interval,
        start_date = _get_start_date_by_wildcards,
        end_date = _get_end_date_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/estimate_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell: """python3 scripts/frequencies.py {input.tree} {input.metadata} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--pivot-frequency {params.pivot_frequency} \
--start-date {params.start_date} \
--end-date {wildcards.timepoint} \
--include-internal-nodes &> {log}"""


rule estimate_diffusion_frequencies:
    message:
        """
        Estimating diffusion frequencies for {input.tree}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        frequencies = BUILD_TIMEPOINT_PATH + "diffusion_frequencies.json"
    params:
        pivot_frequency = _get_pivot_interval,
        stiffness = config["frequencies"]["stiffness"],
        inertia = config["frequencies"]["inertia"],
        min_date = _get_min_date_for_diffusion_frequencies_by_wildcards,
        max_date = _get_max_date_for_augur_frequencies_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_diffusion_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/estimate_diffusion_frequencies_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell: """augur frequencies \
        --method diffusion \
        --tree {input.tree} \
        --metadata {input.metadata} \
        --output {output} \
        --include-internal-nodes \
        --stiffness {params.stiffness} \
        --inertia {params.inertia} \
        --pivot-interval {params.pivot_frequency} \
        --min-date {params.min_date} \
        --max-date {params.max_date} &> {log}"""


rule convert_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_frequencies.output.frequencies
    output:
        table = BUILD_TIMEPOINT_PATH + "frequencies.tsv"
    params:
        method = "kde"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --method {params.method} \
            --output {output} \
            --annotations timepoint={wildcards.timepoint}
        """


rule convert_diffusion_frequencies_to_table:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_diffusion_frequencies.output.frequencies
    output:
        table = BUILD_TIMEPOINT_PATH + "diffusion_frequencies.tsv"
    params:
        method = "diffusion"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/frequencies_to_table.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --method {params.method} \
            --output {output} \
            --annotations timepoint={wildcards.timepoint}
        """


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations for {wildcards}"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment
    output:
        node_data = BUILD_TIMEPOINT_PATH + "nt_muts.json"
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
        reference = "config/reference_h3n2_ha.gb"
    output:
        node_data = BUILD_TIMEPOINT_PATH + "aa_muts.json"
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
        aa_alignment = BUILD_TIMEPOINT_PATH + "aa-seq_{gene}.fasta"
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
        translations = translations(segment="ha", path=BUILD_TIMEPOINT_PATH)
    output:
        translations = BUILD_TIMEPOINT_PATH + "aa_seq.json"
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


rule clades_by_haplotype:
    input:
        tree = rules.refine.output.tree,
        translations = translations(segment="ha", path=BUILD_TIMEPOINT_PATH)
    output:
        clades = BUILD_TIMEPOINT_PATH + "clades.json",
        tip_clade_table = BUILD_TIMEPOINT_PATH + "tips_to_clades.tsv"
    params:
        gene_names = gene_names(segment="ha"),
    conda: "../envs/anaconda.python3.yaml"
    log: "logs/find_clades_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        python3 scripts/nonoverlapping_clades.py \
            --tree {input.tree} \
            --translations {input.translations} \
            --gene-names {params.gene_names} \
            --annotations timepoint={wildcards.timepoint} \
            --output {output.clades} \
            --output-tip-clade-table {output.tip_clade_table} &> {log}
        """


rule delta_frequency:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_diffusion_frequencies.output.frequencies
    output:
        delta_frequency = BUILD_TIMEPOINT_PATH + "delta_frequency.json"
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
            --delta-pivots {params.delta_pivots} \
            --output {output.delta_frequency} &> {log}
        """


rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards
    output:
        node_data = BUILD_TIMEPOINT_PATH + "traits.json",
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
        alignments = translations(segment="ha", path=BUILD_TIMEPOINT_PATH),
        # TODO: define distance maps in build configs
        distance_maps = _get_distance_maps_for_simulations,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names(segment="ha"),
        comparisons = _get_distance_comparisons_for_simulations,
        attribute_names = _get_distance_attributes_for_simulations,
        earliest_date = _get_distance_earliest_date_by_wildcards,
        latest_date = _get_distance_latest_date_by_wildcards
    output:
        distances = BUILD_TIMEPOINT_PATH + "distances.json",
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


rule pairwise_distances:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_frequencies.output.frequencies,
        alignments = translations(segment="ha", path=BUILD_TIMEPOINT_PATH),
        distance_maps = _get_pairwise_distance_maps_for_simulations,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names(segment="ha"),
        attribute_names = _get_pairwise_distance_attributes_for_simulations,
        years_back_to_compare = config["max_years_for_distances"]
    output:
        distances = BUILD_TIMEPOINT_PATH + "pairwise_distances.json",
    benchmark: "benchmarks/pairwise_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/pairwise_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
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


rule cross_immunities:
    input:
        frequencies = rules.estimate_frequencies.output.frequencies,
        distances = rules.pairwise_distances.output.distances,
        date_annotations = rules.refine.output.node_data
    params:
        distance_attributes = _get_cross_immunity_distance_attributes_for_simulations,
        immunity_attributes = _get_cross_immunity_attributes_for_simulations,
        decay_factors = _get_cross_immunity_decay_factors_for_simulations,
        years_to_wane = config["max_years_for_distances"]
    output:
        cross_immunities = BUILD_TIMEPOINT_PATH + "cross_immunity.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/cross_immunity.py \
            --frequencies {input.frequencies} \
            --distances {input.distances} \
            --date-annotations {input.date_annotations} \
            --distance-attributes {params.distance_attributes} \
            --immunity-attributes {params.immunity_attributes} \
            --decay-factors {params.decay_factors} \
            --years-to-wane {params.years_to_wane} \
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
        lbi = BUILD_TIMEPOINT_PATH + "lbi.json"
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


rule unnormalized_lbi:
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = config["lbi"]["tau"],
        window = config["lbi"]["window"],
        names = "unnormalized_lbi"
    output:
        lbi = BUILD_TIMEPOINT_PATH + "unnormalized_lbi.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window} \
            --no-normalization
        """


rule filter_translations_by_date:
    input:
        alignments = rules.reconstruct_translations.output.aa_alignment,
        branch_lengths = rules.refine.output.node_data
    output:
        alignments = BUILD_TIMEPOINT_PATH + "filtered-aa-seq_{gene}.fasta"
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
        titers = _get_titers_by_wildcards,
        alignments = titer_sequence_translations
    params:
        genes = gene_names
    output:
        titers_model = BUILD_TIMEPOINT_PATH + "titers-sub-model.json",
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
        titers = _get_titers_by_wildcards,
        tree = rules.refine.output.tree
    output:
        titers_model = BUILD_TIMEPOINT_PATH + "titers-tree-model.json",
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
        distance_map = BUILD_TIMEPOINT_PATH + "titer_substitution_distance_map.json"
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
        # TODO: move these params to builds in config file
        genes = gene_names,
        comparisons = "root ancestor",
        attribute_names = "cTiterSub cTiterSub_star",
        earliest_date = _get_distance_earliest_date_by_wildcards,
        latest_date = _get_distance_latest_date_by_wildcards
    output:
        distances = BUILD_TIMEPOINT_PATH + "titer_substitution_distances.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --compare-to {params.comparisons} \
            --attribute-name {params.attribute_names} \
            --map {input.distance_maps} {input.distance_maps} \
            --date-annotations {input.date_annotations} \
            --earliest-date {params.earliest_date} \
            --latest-date {params.latest_date} \
            --output {output}
        """


rule pairwise_titer_distances:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.estimate_frequencies.output.frequencies,
        alignments = translations(segment="ha", path=BUILD_TIMEPOINT_PATH),
        distance_maps = rules.convert_titer_model_to_distance_map.output.distance_map,
        date_annotations = rules.refine.output.node_data
    params:
        genes = gene_names(segment="ha"),
        attribute_names = "cTiterSub_pairwise",
        years_back_to_compare = config["max_years_for_distances"]
    output:
        distances = BUILD_TIMEPOINT_PATH + "pairwise_titer_distances.json",
    benchmark: "benchmarks/pairwise_titer_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/pairwise_titer_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
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


rule titer_cross_immunities:
    input:
        frequencies = rules.estimate_frequencies.output.frequencies,
        distances = rules.pairwise_titer_distances.output.distances,
        date_annotations = rules.refine.output.node_data
    params:
        distance_attributes = "cTiterSub_pairwise",
        immunity_attributes = "cTiterSub_x",
        decay_factors = "14.0",
        years_to_wane = config["max_years_for_distances"]
    output:
        cross_immunities = BUILD_TIMEPOINT_PATH + "titer_substitution_cross_immunity.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/cross_immunity.py \
            --frequencies {input.frequencies} \
            --distances {input.distances} \
            --date-annotations {input.date_annotations} \
            --distance-attributes {params.distance_attributes} \
            --immunity-attributes {params.immunity_attributes} \
            --decay-factors {params.decay_factors} \
            --years-to-wane {params.years_to_wane} \
            --output {output}
        """


rule normalize_fitness:
    input:
        metadata = _get_metadata_by_wildcards,
        frequencies = rules.convert_frequencies_to_table.output.table
    output:
        fitness = BUILD_TIMEPOINT_PATH + "normalized_fitness.json"
    params:
        preferred_frequency_method = config["frequencies"]["preferred_method"]
    shell:
        """
        python3 scripts/normalize_fitness.py \
            --metadata {input.metadata} \
            --frequencies-table {input.frequencies} \
            --frequency-method {params.preferred_frequency_method} \
            --output {output.fitness}
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
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        weights = "data/region_weights.json"
    output:
        frequencies = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_tip-frequencies.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=_get_pivot_interval,
        min_date=_get_min_date_for_augur_frequencies_by_wildcards,
        max_date=_get_max_date_for_augur_frequencies_by_wildcards
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
        rules.clades_by_haplotype.output.clades,
        rules.delta_frequency.output.delta_frequency,
        rules.distances.output.distances,
        rules.cross_immunities.output.cross_immunities,
        rules.lbi.output.lbi,
        rules.unnormalized_lbi.output.lbi
    ]

    # Define node data that only make sense for natural populations
    # such as titer models.
    if wildcards.type == "natural":
        inputs.extend([
            rules.traits.output.node_data,
            rules.titers_tree.output.titers_model,
            rules.titer_distances.output.distances,
            rules.titer_cross_immunities.output.cross_immunities
        ])
    elif wildcards.type == "simulated":
        inputs.extend([
            rules.normalize_fitness.output.fitness
        ])

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs


rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        auspice_config = "config/auspice_config.json",
        node_data = _get_node_data_for_export,
        colors = "config/colors.tsv"
    output:
        auspice_tree = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_tree.json",
        auspice_metadata = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_meta.json"
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


rule convert_node_data_to_table:
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        node_data = _get_node_data_for_export
    output:
        table = BUILD_TIMEPOINT_PATH + "node_data.tsv"
    params:
        excluded_fields_arg = _get_excluded_fields_arg,
        lineage = _get_lineage,
        segment = _get_segment
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
                          lineage={params.lineage} \
                          segment={params.segment}
        """


rule merge_node_data_and_frequencies:
    input:
        node_data = rules.convert_node_data_to_table.output.table,
        kde_frequencies = rules.convert_frequencies_to_table.output.table,
        diffusion_frequencies = rules.convert_diffusion_frequencies_to_table.output.table
    output:
        table = BUILD_TIMEPOINT_PATH + "tip_attributes.tsv"
    run:
        node_data = pd.read_table(input.node_data)
        kde_frequencies = pd.read_table(input.kde_frequencies)
        diffusion_frequencies = pd.read_table(input.diffusion_frequencies)
        df = node_data.merge(
            kde_frequencies,
            how="inner",
            on=["strain", "timepoint", "is_terminal"]
        ).merge(
            diffusion_frequencies,
            how="inner",
            on=["strain", "timepoint", "is_terminal"]
        )

        df.to_csv(output.table, sep="\t", index=False, header=True)


rule collect_tip_attributes:
    input:
        _get_tip_attributes_by_wildcards
    output:
        attributes = BUILD_PATH + "tip_attributes.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/collect_tables.py \
            --tables {input} \
            --output {output.attributes}
        """


rule annotate_naive_tip_attribute:
    input:
        attributes = rules.collect_tip_attributes.output.attributes
    output:
        attributes = BUILD_PATH + "tip_attributes_with_naive_predictor.tsv",
    params:
        preferred_frequency_method = config["frequencies"]["preferred_method"]
    run:
        # Annotate a predictor for a naive model with no growth.
        df = pd.read_csv(input.attributes, sep="\t")
        df["naive"] = 0.0
        df["frequency"] = df["%s_frequency" % params.preferred_frequency_method]
        df = df[df["frequency"] > 0.0].copy()
        df.to_csv(output.attributes, sep="\t", index=False)


rule annotate_observed_offspring:
    input:
        tree = _get_final_tree_for_wildcards,
        attributes = rules.annotate_naive_tip_attribute.output.attributes
    output:
        attributes = BUILD_PATH + "tip_attributes_with_offspring.tsv",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/annotate_offspring.py \
            --tree {input.tree} \
            --tip-attributes {input.attributes} \
            --output {output.attributes}
        """


rule target_distances:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes
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


rule annotate_weighted_distances_for_tip_attributes:
    input:
        attributes = rules.annotate_observed_offspring.output.attributes,
        distances = rules.target_distances.output.distances
    output:
        attributes = BUILD_PATH + "tip_attributes_with_weighted_distances.tsv"
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


rule fit_models_by_distances:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        distances = rules.target_distances.output.distances
    output:
        model = BUILD_PATH + "models_by_distances/{predictors}.json",
        errors = BUILD_PATH + "models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "models_by_distances_coefficients/{predictors}.tsv"
    params:
        predictors = _get_predictor_list,
        delta_months = config["fitness_model"]["delta_months"],
        training_window = config["fitness_model"]["training_window"],
        cost_function = config["fitness_model"]["distance_cost_function"],
        l1_lambda = config["fitness_model"]["l1_lambda"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/fitness_model_distances_" + BUILD_LOG_STEM + "_{predictors}.txt"
    log: "logs/fitness_model_distances_" + BUILD_LOG_STEM + "_{predictors}.txt"
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
            --errors-by-timepoint {output.errors} \
            --coefficients-by-timepoint {output.coefficients} \
            --include-scores \
            --output {output.model} &> {log}
        """


rule annotate_distance_models:
    input:
        errors = rules.fit_models_by_distances.output.errors,
        coefficients = rules.fit_models_by_distances.output.coefficients
    output:
        errors = BUILD_PATH + "annotated_models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "annotated_models_by_distances_coefficients/{predictors}.tsv"
    run:
        errors = pd.read_csv(input.errors, sep="\t")
        errors["type"] = wildcards.type
        errors["sample"] = wildcards.sample
        errors.to_csv(output.errors, sep="\t", header=True, index=False)

        coefficients = pd.read_csv(input.coefficients, sep="\t")
        coefficients["type"] = wildcards.type
        coefficients["sample"] = wildcards.sample
        coefficients.to_csv(output.coefficients, sep="\t", header=True, index=False)


rule collect_annotated_tip_clade_tables:
    input:
        _get_tip_clades_by_wildcards
    output:
        tip_clade_table = BUILD_PATH + "tips_to_clades.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/collect_tables.py \
            --tables {input} \
            --output {output.tip_clade_table}
        """


rule select_clades:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        tips_to_clades = rules.collect_annotated_tip_clade_tables.output.tip_clade_table
    output:
        clades = BUILD_PATH + "final_clade_frequencies.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months"]
    conda: "../envs/anaconda.python3.yaml"
    log: "logs/select_clades_" + BUILD_LOG_STEM + ".txt"
    shell:
        """
        python3 scripts/select_clades.py \
            --tip-attributes {input.attributes} \
            --tips-to-clades {input.tips_to_clades} \
            --delta-months {params.delta_months} \
            --output {output} &> {log}
        """


rule fit_models_by_clades:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes,
        final_clade_frequencies = rules.select_clades.output.clades
    output:
        model = BUILD_PATH + "models_by_clades/{predictors}.json"
    params:
        predictors = _get_predictor_list,
        delta_months = config["fitness_model"]["delta_months"],
        training_window = config["fitness_model"]["training_window"],
        cost_function = config["fitness_model"]["clade_cost_function"],
        l1_lambda = config["fitness_model"]["l1_lambda"],
        pseudocount = config["fitness_model"]["pseudocount"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/fitness_model_clades_" + BUILD_LOG_STEM + "_{predictors}.txt"
    log: "logs/fitness_model_clades_" + BUILD_LOG_STEM + "_{predictors}.txt"
    shell:
        """
        python3 src/fit_model.py \
            --tip-attributes {input.attributes} \
            --final-clade-frequencies {input.final_clade_frequencies} \
            --training-window {params.training_window} \
            --delta-months {params.delta_months} \
            --predictors {params.predictors} \
            --cost-function {params.cost_function} \
            --l1-lambda {params.l1_lambda} \
            --pseudocount {params.pseudocount} \
            --target clades \
            --output {output} &> {log}
        """


# rule export_with_fitness:
#     input:
#         tree = rules.refine.output.tree,
#         metadata = _get_metadata_by_wildcards,
#         auspice_config = "config/auspice_config.json",
#         node_data = _get_node_data_for_export,
#         fitness = BUILD_TIMEPOINT_PATH + "fitness.json",
#         colors = "config/colors.tsv"
#     output:
#         auspice_tree = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_tree.json",
#         auspice_metadata = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_meta.json"
#     params:
#         panels = "tree entropy frequencies"
#     conda: "../envs/anaconda.python3.yaml"
#     shell:
#         """
#         augur export \
#             --tree {input.tree} \
#             --metadata {input.metadata} \
#             --node-data {input.node_data} \
#             --colors {input.colors} \
#             --auspice-config {input.auspice_config} \
#             --output-tree {output.auspice_tree} \
#             --output-meta {output.auspice_metadata} \
#             --panels {params.panels} \
#             --minify-json
#         """


# rule plot_tree:
#     input:
#         auspice_tree = rules.export.output.auspice_tree
#     output:
#         tree = "results/figures/trees/flu_" + BUILD_SEGMENT_LOG_STEM + "_tree.pdf"
#     conda: "../envs/anaconda.python3.yaml"
#     benchmark: "benchmarks/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
#     log: "logs/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
#     params:
#         start = _get_start_date_by_wildcards,
#         end = _get_end_date_by_wildcards
#     shell:
#         """
#         python3 scripts/plot_tree.py \
#             {input} \
#             {output} \
#             --start-date {params.start} \
#             --end-date {params.end} &> {log}
#         """

# # rule aggregate_tree_plots_simulated:
# #     input: expand(rules.plot_tree_simulated.output.tree, percentage=PERCENTAGE, start=START_DATE_SIMULATIONS, end=END_DATE_SIMULATIONS, timepoint=TIMEPOINTS_SIMULATIONS)
# #     output:
# #         trees="results/figures/trees_simulated.pdf"
# #     shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"
