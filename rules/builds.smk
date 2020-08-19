"""Rules for generating simulated HA sequences for validation of forecasting models.
"""
BUILD_PATH = "results/builds/{type}/{sample}/"
BUILD_LOG_STEM = "{type}_{sample}"
BUILD_TIMEPOINT_PATH = BUILD_PATH + "timepoints/{timepoint}/"
BUILD_SEGMENT_LOG_STEM = "{type}_{sample}_{timepoint}"


rule get_strains_by_timepoint:
    input:
        metadata = _get_metadata_by_wildcards
    output:
        strains = BUILD_TIMEPOINT_PATH + "strains.txt"
    params:
        years_back = _get_years_back_to_build_trees,
        reference_strains = _get_required_strains_argument
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/partition_strains_by_timepoint.py \
            {input.metadata} \
            {wildcards.timepoint} \
            {output} \
            --years-back {params.years_back} \
            {params.reference_strains}
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
        metadata = _get_metadata_by_wildcards
    output:
        frequencies = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_original-tip-frequencies.json"
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
        reference = _get_reference
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
    return config["cross_immunity"]["h3n2"]["ha"][wildcards.type]["distance_attributes"]


def _get_cross_immunity_attributes_for_simulations(wildcards):
    return config["cross_immunity"]["h3n2"]["ha"][wildcards.type]["immunity_attributes"]


def _get_cross_immunity_decay_factors_for_simulations(wildcards):
    return config["cross_immunity"]["h3n2"]["ha"][wildcards.type]["decay_factors"]


rule cross_immunities:
    input:
        frequencies = rules.tip_frequencies.output.frequencies,
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


# Only attempt to regenerate titer model files if the user has access to the raw titer data.
# Otherwise, the files stored in version control will be used.
if "RETHINK_HOST" in os.environ and "RETHINK_AUTH_KEY" in os.environ:
    rule titers_sub:
        input:
            titers = _get_titers_by_wildcards,
            alignments = translations,
            tree = rules.refine.output.tree
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
                --tree {input.tree} \
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


    rule fra_titers_tree:
        input:
            titers = _get_fra_titers_by_wildcards,
            tree = rules.refine.output.tree
        output:
            titers_model = BUILD_TIMEPOINT_PATH + "fra-titers-tree-model.json",
        conda: "../envs/anaconda.python3.yaml"
        benchmark: "benchmarks/fra_titers_tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
        log: "logs/fra_titers_tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
        shell:
            """
            augur titers tree \
                --titers {input.titers} \
                --tree {input.tree} \
                --allow-empty-model \
                --output {output.titers_model} &> {log}
            """


# This is silly, but augur titers outputs a fixed pair of key names and trying
# to merge two different titer models into an auspice JSON will cause a
# collision.
rule rename_fields_in_fra_titers_tree:
    input:
        titers_model = BUILD_TIMEPOINT_PATH + "fra-titers-tree-model.json"
    output:
        titers_model = BUILD_TIMEPOINT_PATH + "renamed-fra-titers-tree-model.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/rename_fields_in_fra_titer_models.py \
            --titers-model {input.titers_model} \
            --output {output.titers_model}
        """


rule convert_titer_model_to_distance_map:
    input:
        model = BUILD_TIMEPOINT_PATH + "titers-sub-model.json"
    output:
        distance_map = BUILD_TIMEPOINT_PATH + "titer_substitution_distance_map.json"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/titer_model_to_distance_map.py \
            --model {input.model} \
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


rule pairwise_titer_tree_distances:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.tip_frequencies.output.frequencies,
        model = BUILD_TIMEPOINT_PATH + "titers-tree-model.json",
        date_annotations = rules.refine.output.node_data
    params:
        attribute_names = "cTiter_pairwise",
        months_back_for_current_samples = config["months_back_for_current_samples"],
        years_back_to_compare = config["max_years_for_distances"]
    output:
        distances = BUILD_TIMEPOINT_PATH + "pairwise_titer_tree_distances.json",
    benchmark: "benchmarks/pairwise_titer_tree_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/pairwise_titer_tree_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/pairwise_titer_tree_distances.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --attribute-name {params.attribute_names} \
            --date-annotations {input.date_annotations} \
            --months-back-for-current-samples {params.months_back_for_current_samples} \
            --years-back-to-compare {params.years_back_to_compare} \
            --output {output} &> {log}
        """


rule pairwise_fra_titer_tree_distances:
    input:
        tree = rules.refine.output.tree,
        frequencies = rules.tip_frequencies.output.frequencies,
        model = rules.rename_fields_in_fra_titers_tree.output.titers_model,
        date_annotations = rules.refine.output.node_data
    params:
        attribute_names = "fra_cTiter_pairwise",
        model_attribute_name = "fra_dTiter",
        months_back_for_current_samples = config["months_back_for_current_samples"],
        years_back_to_compare = config["max_years_for_distances"],
    output:
        distances = BUILD_TIMEPOINT_PATH + "pairwise_fra_titer_tree_distances.json",
    benchmark: "benchmarks/pairwise_fra_titer_tree_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/pairwise_fra_titer_tree_distances_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/pairwise_titer_tree_distances.py \
            --tree {input.tree} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --model-attribute-name {params.model_attribute_name} \
            --attribute-name {params.attribute_names} \
            --date-annotations {input.date_annotations} \
            --months-back-for-current-samples {params.months_back_for_current_samples} \
            --years-back-to-compare {params.years_back_to_compare} \
            --output {output} &> {log}
        """


rule titer_cross_immunities:
    input:
        frequencies = rules.tip_frequencies.output.frequencies,
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
            --output {output}
        """


rule titer_tree_cross_immunities:
    input:
        frequencies = rules.tip_frequencies.output.frequencies,
        distances = rules.pairwise_titer_tree_distances.output.distances,
        date_annotations = rules.refine.output.node_data
    params:
        distance_attributes = "cTiter_pairwise",
        immunity_attributes = "cTiter_x",
        decay_factors = "14.0",
        years_to_wane = config["max_years_for_distances"]
    output:
        cross_immunities = BUILD_TIMEPOINT_PATH + "titer_tree_cross_immunity.json",
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
            --output {output}
        """


rule fra_titer_tree_cross_immunities:
    input:
        frequencies = rules.tip_frequencies.output.frequencies,
        distances = rules.pairwise_fra_titer_tree_distances.output.distances,
        date_annotations = rules.refine.output.node_data
    params:
        distance_attributes = "fra_cTiter_pairwise",
        immunity_attributes = "fra_cTiter_x",
        decay_factors = "14.0",
        years_to_wane = config["max_years_for_distances"]
    output:
        cross_immunities = BUILD_TIMEPOINT_PATH + "fra_titer_tree_cross_immunity.json",
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
        rules.lbi.output.lbi
    ]

    # Define node data that only make sense for natural populations
    # such as titer models.
    if wildcards.type == "natural":
        inputs.extend([
            rules.traits.output.node_data,
            BUILD_TIMEPOINT_PATH + "titers-tree-model.json",
            BUILD_TIMEPOINT_PATH + "titers-sub-model.json",
            rules.titer_cross_immunities.output.cross_immunities,
            rules.titer_tree_cross_immunities.output.cross_immunities
        ])

        build = config["builds"][wildcards.type][wildcards.sample]
        if "fra_titers" in build:
           inputs.append(rules.rename_fields_in_fra_titers_tree.output.titers_model)
           inputs.append(rules.fra_titer_tree_cross_immunities.output.cross_immunities)

    elif wildcards.type == "simulated":
        inputs.extend([
            rules.normalize_fitness.output.fitness
        ])

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs


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
    params:
        preferred_frequency_method = config["frequencies"]["preferred_method"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/merge_node_data_and_frequencies.py \
            --node-data {input.node_data} \
            --kde-frequencies {input.kde_frequencies} \
            --diffusion-frequencies {input.diffusion_frequencies} \
            --preferred-frequency-method {params.preferred_frequency_method} \
            --output {output.table}
        """


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
        attributes = BUILD_PATH + "tip_attributes_with_naive_predictor.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/annotate_naive_tip_attribute.py \
            --tip-attributes {input.attributes} \
            --output {output.attributes}
        """


rule target_distances:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes
    output:
        distances = BUILD_PATH + "target_distances.tsv",
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        sequence_attribute_name = "aa_sequence"
    benchmark: "benchmarks/target_distances_" + BUILD_LOG_STEM + ".txt"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --sequence-attribute-name {params.sequence_attribute_name} \
            --output {output}
        """


rule annotate_weighted_distances_for_tip_attributes:
    input:
        attributes = rules.annotate_naive_tip_attribute.output.attributes,
        distances = rules.target_distances.output.distances
    output:
        attributes = BUILD_PATH + "tip_attributes_with_weighted_distances.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"]
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
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        training_window = _get_fitness_model_training_window,
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


rule extract_minimal_models_by_distances:
    input:
        model = rules.fit_models_by_distances.output.model
    output:
        model = BUILD_PATH + "minimal_models_by_distances/{predictors}.json",
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/extract_minimal_models_by_distances.py \
            --model {input.model} \
            --output {output.model}
        """


rule annotate_distance_models:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        model = rules.fit_models_by_distances.output.model,
        errors = rules.fit_models_by_distances.output.errors,
        coefficients = rules.fit_models_by_distances.output.coefficients
    output:
        errors = BUILD_PATH + "annotated_models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "annotated_models_by_distances_coefficients/{predictors}.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        error_type = "validation"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/annotate_model_tables.py \
            --tip-attributes {input.attributes} \
            --model {input.model} \
            --errors-by-timepoint {input.errors} \
            --coefficients-by-timepoint {input.coefficients} \
            --annotated-errors-by-timepoint {output.errors} \
            --annotated-coefficients-by-timepoint {output.coefficients} \
            --delta-months {params.delta_months} \
            --annotations type="{wildcards.type}" sample="{wildcards.sample}" error_type="{params.error_type}"
        """


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
        delta_months = config["fitness_model"]["delta_months_to_fit"]
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
        delta_months = config["fitness_model"]["delta_months_to_fit"],
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


rule plot_tree:
    input:
        tree = rules.refine.output.tree
    output:
        tree = BUILD_TIMEPOINT_PATH + "tree.pdf"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        python3 scripts/plot_tree.py {input} {output} &> {log}
        """


rule aggregate_tree_plots:
    input: _get_tree_plots_by_wildcards
    output:
        trees="results/figures/trees_" + BUILD_LOG_STEM + ".pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"


rule target_distances_by_timepoint:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table
    output:
        distances = BUILD_TIMEPOINT_PATH + "target_distances.tsv",
    params:
        delta_months = _get_delta_months_to_forecast
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/calculate_target_distances.py \
            --tip-attributes {input.attributes} \
            --delta-months {params.delta_months} \
            --output {output}
        """


rule forecast_tips:
    input:
        attributes = rules.merge_node_data_and_frequencies.output.table,
        distances = rules.target_distances_by_timepoint.output.distances,
        frequencies = rules.tip_frequencies.output.frequencies,
        model = lambda wildcards: config["builds"][wildcards.type][wildcards.sample]["best_predictor"]
    output:
        node_data = BUILD_TIMEPOINT_PATH + "forecasts.json",
        frequencies = "results/auspice/flu_" + BUILD_SEGMENT_LOG_STEM + "_tip-frequencies.json"
    params:
        delta_months = _get_delta_months_to_forecast
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/forecast_model.py \
            --tip-attributes {input.attributes} \
            --distances {input.distances} \
            --frequencies {input.frequencies} \
            --model {input.model} \
            --delta-months {params.delta_months} \
            --output-node-data {output.node_data} \
            --output-frequencies {output.frequencies}
        """


rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = _get_metadata_by_wildcards,
        auspice_config = "config/auspice_config.json",
        node_data = _get_node_data_for_export,
        forecasts = rules.forecast_tips.output.node_data,
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
            --node-data {input.node_data} {input.forecasts} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_metadata} \
            --panels {params.panels} \
            --minify-json
        """


rule forecast_all_tips:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        distances = rules.target_distances.output.distances,
        model = _get_model_from_validation
    output:
        table = BUILD_PATH + "forecasts_{predictors}.tsv",
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"]
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 src/forecast_model.py \
            --tip-attributes {input.attributes} \
            --distances {input.distances} \
            --model {input.model} \
            --delta-months {params.delta_months} \
            --output-table {output.table}
        """


rule test_distance_models:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        distances = rules.target_distances.output.distances,
        model = _get_model_to_test_by_wildcards
    output:
        model = BUILD_PATH + "test_models_by_distances/{predictors}.json",
        errors = BUILD_PATH + "test_models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "test_models_by_distances_coefficients/{predictors}.tsv"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/test_fitness_model_distances_" + BUILD_LOG_STEM + "_{predictors}.txt"
    log: "logs/test_fitness_model_distances_" + BUILD_LOG_STEM + "_{predictors}.txt"
    shell:
        """
        python3 src/fit_model.py \
            --tip-attributes {input.attributes} \
            --target distances \
            --distances {input.distances} \
            --fixed-model {input.model} \
            --errors-by-timepoint {output.errors} \
            --coefficients-by-timepoint {output.coefficients} \
            --include-scores \
            --output {output.model} &> {log}
        """


rule annotate_test_distance_models:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        model = rules.test_distance_models.output.model,
        errors = rules.test_distance_models.output.errors,
        coefficients = rules.test_distance_models.output.coefficients
    output:
        errors = BUILD_PATH + "annotated_test_models_by_distances_errors/{predictors}.tsv",
        coefficients = BUILD_PATH + "annotated_test_models_by_distances_coefficients/{predictors}.tsv"
    params:
        delta_months = config["fitness_model"]["delta_months_to_fit"],
        sample = _get_validation_sample_by_wildcards
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/annotate_model_tables.py \
            --tip-attributes {input.attributes} \
            --model {input.model} \
            --errors-by-timepoint {input.errors} \
            --coefficients-by-timepoint {input.coefficients} \
            --annotated-errors-by-timepoint {output.errors} \
            --annotated-coefficients-by-timepoint {output.coefficients} \
            --delta-months {params.delta_months} \
            --annotations type="{wildcards.type}" sample="{params.sample}" error_type="test"
        """


def _get_tips_to_clades_for_full_tree_by_wildcards(wildcards):
    full_tree_sample = _get_full_tree_sample_by_wildcards(wildcards)
    path = rules.collect_annotated_tip_clade_tables.output.tip_clade_table.format(
        type=wildcards.type,
        sample=full_tree_sample
    )
    return path


rule plot_validation_figure:
    input:
        attributes = rules.annotate_weighted_distances_for_tip_attributes.output.attributes,
        tips_to_clades = _get_tips_to_clades_for_full_tree_by_wildcards,
        forecasts = rules.forecast_all_tips.output.table,
        model_errors = _get_model_validation_errors
    output:
        figure = "manuscript/figures/validation_figure_{type}-{sample}-{predictors}.pdf",
        clades = "manuscript/figures/validation_figure_clades_{type}-{sample}-{predictors}.tsv",
        ranks = "manuscript/figures/validation_figure_ranks_{type}-{sample}-{predictors}.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/plot_validation_figure_by_population.py \
            --tip-attributes {input.attributes} \
            --tips-to-clades {input.tips_to_clades} \
            --forecasts {input.forecasts} \
            --model-errors {input.model_errors} \
            --population {wildcards.type} \
            --sample {wildcards.sample} \
            --predictors {wildcards.predictors} \
            --output {output.figure} \
            --output-clades-table {output.clades} \
            --output-ranks-table {output.ranks}
        """
