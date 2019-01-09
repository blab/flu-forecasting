"""
Rules for quality control visualizations including metadata, trees, and frequencies.
"""

#
# Visualize frequencies
#

rule plot_frequencies:
    input: rules.estimate_frequencies.output.frequencies
    output: "results/figures/frequencies/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python3 scripts/plot_frequency_trajectories.py {input} {output}"

rule aggregate_frequency_plots:
    input: expand("results/figures/frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "results/figures/frequencies.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

#
# Visualize sequence metadata
#

rule plot_sequences_by_date:
    input: rules.export.output.auspice_tree
    output: "results/figures/sequence_distributions/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python3 scripts/plot_sequence_distribution.py {input} {output}"

rule aggregate_sequence_distribution_plots:
    input: expand("results/figures/sequence_distributions/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "results/figures/sequence_distributions.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

#
# Visualize trees
#

rule plot_tree:
    input: rules.export.output.auspice_tree
    output: "results/figures/trees/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}_tree.pdf"
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/plot_tree_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/plot_tree_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.log"
    shell: """python3 scripts/plot_tree.py {input} {output} &> {log}"""

rule aggregate_tree_plots:
    input: expand("results/figures/trees/flu_h3n2_{segment}_{year_range}y_{viruses}v_{sample}_tree.pdf", segment=SEGMENTS, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output:
        trees="results/figures/trees.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"
