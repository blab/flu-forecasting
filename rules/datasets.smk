"""
Rules to build simulated datasets.
"""
from augur.frequency_estimators import float_to_datestring
import Bio.SeqIO
import os
import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

DATA_SIMULATED_ROOT_PATH = "data/{type}/{sample}/"


def _get_viruses_per_month(wildcards):
    return config["datasets"][wildcards.type][wildcards.sample]["viruses_per_month"]

def _get_simulation_seed(wildcards):
    return config["datasets"][wildcards.type][wildcards.sample]["seed"]


rule run_simulation:
    input:
        simulation_config = "data/{type}/{sample}/influenza_h3n2_ha.xml"
    output:
        sequences = "data/{type}/{sample}/simulated_HA_sequences.fasta"
    params:
        seed = _get_simulation_seed
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        cd data/{wildcards.type}/{wildcards.sample} && java -jar {SNAKEMAKE_DIR}/dist/santa-sim/dist/santa.jar -seed={params.seed} {SNAKEMAKE_DIR}/{input.simulation_config}
        """


rule parse_simulated_sequences:
    input:
        sequences = rules.run_simulation.output.sequences
    output:
        sequences = DATA_SIMULATED_ROOT_PATH + "sequences.fasta",
        metadata = DATA_SIMULATED_ROOT_PATH + "metadata.tsv"
    params:
        fasta_fields = "strain generation fitness"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """


rule standardize_simulated_sequence_dates:
    input:
        metadata = rules.parse_simulated_sequences.output.metadata
    output:
        metadata = DATA_SIMULATED_ROOT_PATH + "corrected_metadata.tsv"
    run:
        df = pd.read_csv(input.metadata, sep="\t")
        df["num_date"] = 2000.0 + (df["generation"] / 100.0)
        df["date"] = df["num_date"].apply(float_to_datestring)
        df["year"]  = pd.to_datetime(df["date"]).dt.year
        df["month"]  = pd.to_datetime(df["date"]).dt.month

        df[df["fitness"] > 0].to_csv(output.metadata, header=True, index=False, sep="\t")


rule filter_simulated:
    input:
        sequences = rules.parse_simulated_sequences.output.sequences,
        metadata = rules.standardize_simulated_sequence_dates.output.metadata
    output:
        sequences = DATA_SIMULATED_ROOT_PATH + "filtered_sequences.fasta"
    params:
        # Skip the first 1,000 generations (or 1000 / 100 years) for simulation burn-in.
        min_date = 2010.0,
        group_by = "year month",
        viruses_per_month = _get_viruses_per_month
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/filter_{type}_{sample}.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-date {params.min_date} \
            --group-by {params.group_by} \
            --sequences-per-group {params.viruses_per_month} \
            --output {output}
        """


rule filter_metadata_simulated:
    input:
        sequences = rules.filter_simulated.output.sequences,
        metadata = rules.standardize_simulated_sequence_dates.output.metadata,
    output:
        metadata = DATA_SIMULATED_ROOT_PATH + "filtered_metadata.tsv"
    run:
        # Get a list of all samples that passed the sequence filtering step.
        sequences = Bio.SeqIO.parse(input.sequences, "fasta")
        sample_ids = [sequence.id for sequence in sequences]

        # Load all metadata.
        metadata = pd.read_csv(input.metadata, sep="\t")
        filtered_metadata = metadata[metadata["strain"].isin(sample_ids)].copy()

        # Save only the metadata records that have entries in the filtered sequences.
        filtered_metadata.to_csv(output.metadata, sep="\t", header=True, index=False)


rule get_strains_for_simulated_sequences:
    input:
        metadata = rules.filter_metadata_simulated.output.metadata
    output:
        strains = DATA_SIMULATED_ROOT_PATH + "strains.txt"
    run:
        df = pd.read_csv(input.metadata, sep="\t")
        df["strain"].to_csv(output.strains, header=False, index=False)
