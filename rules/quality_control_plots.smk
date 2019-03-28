"""
Rules for quality control visualizations including metadata, trees, and frequencies.
"""

#
# Visualize frequencies
#

rule plot_frequencies:
    input:
        frequencies = rules.estimate_frequencies.output.frequencies
    output:
        frequencies = "results/figures/frequencies/flu_" + BUILD_SEGMENT_LOG_STEM + ".pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python3 scripts/plot_frequency_trajectories.py {input} {output}"

rule aggregate_frequency_plots:
    input: expand(rules.plot_frequencies.output.frequencies, lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, timepoint=TIMEPOINTS, segment=SEGMENTS)
    output: "results/figures/frequencies.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

#
# Visualize sequence metadata
#

rule plot_sequences_by_date:
    input:
        auspice_tree = rules.export.output.auspice_tree
    output:
        sequences = "results/figures/sequence_distributions/flu_" + BUILD_SEGMENT_LOG_STEM + ".pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python3 scripts/plot_sequence_distribution.py {input} {output}"

rule aggregate_sequence_distribution_plots:
    input: expand(rules.plot_sequences_by_date.output.sequences, lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, timepoint=TIMEPOINTS, segment=SEGMENTS)
    output: "results/figures/sequence_distributions.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

#
# Visualize trees
#

rule plot_tree:
    input:
        auspice_tree = rules.export.output.auspice_tree
    output:
        tree = "results/figures/trees/flu_" + BUILD_SEGMENT_LOG_STEM + "_tree.pdf"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".txt"
    log: "logs/plot_tree_" + BUILD_SEGMENT_LOG_STEM + ".log"
    shell:
        """
        python3 scripts/plot_tree.py \
            {input} \
            {output} \
            --start-date {wildcards.start} \
            --end-date {wildcards.end} &> {log}
        """

rule aggregate_tree_plots:
    input: expand(rules.plot_tree.output.tree, lineage=LINEAGES, viruses=VIRUSES, sample=SAMPLES, start=START_DATE, end=END_DATE, timepoint=TIMEPOINTS, segment=SEGMENTS)
    output:
        trees="results/figures/trees.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"
