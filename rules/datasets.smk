"""
Rules to build simulated datasets.
"""
# Set template for root path of any given dataset.
# Datasets are defined by their type (e.g., natural, simulated, etc.) and their sample name.
DATA_SIMULATED_ROOT_PATH = "data/simulated/{sample}/"
DATA_NATURAL_ROOT_PATH = "data/natural/{sample}/"

#
# Rules for simulated datasets
#


rule run_simulation:
    input:
        simulation_config = ancient(DATA_SIMULATED_ROOT_PATH + "influenza_h3n2_ha.xml")
    output:
        sequences = DATA_SIMULATED_ROOT_PATH + "simulated_HA_sequences.fasta"
    params:
        seed = _get_simulation_seed
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        cd data/simulated/{wildcards.sample} && java -jar {SNAKEMAKE_DIR}/dist/santa.jar -seed={params.seed} {SNAKEMAKE_DIR}/{input.simulation_config}
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
    conda: "../envs/anaconda.python3.yaml"
    params:
        start_year = 2000.0,
        generations_per_year = 200.0
    shell:
        """
        python3 scripts/standardize_simulated_sequence_dates.py \
            --metadata {input.metadata} \
            --start-year {params.start_year} \
           --generations-per-year {params.generations_per_year} \
            --output {output.metadata}
        """


rule filter_simulated:
    input:
        sequences = ancient(rules.parse_simulated_sequences.output.sequences),
        metadata = ancient(rules.standardize_simulated_sequence_dates.output.metadata)
    output:
        sequences = DATA_SIMULATED_ROOT_PATH + "filtered_sequences.fasta"
    params:
        # Skip the first 1,000 generations (or 1000 / 100 years) for simulation burn-in.
        min_date = 2010.0,
        group_by = "year month",
        viruses_per_month = _get_viruses_per_month
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/filter_simulated_{sample}.txt"
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
        sequences = ancient(rules.filter_simulated.output.sequences),
        metadata = ancient(rules.standardize_simulated_sequence_dates.output.metadata)
    output:
        metadata = DATA_SIMULATED_ROOT_PATH + "filtered_metadata.tsv"
    conda: "../envs/anaconda.python3.yaml"
    shell:
        """
        python3 scripts/filter_simulated_metadata.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output {output.metadata}
        """


#
# Rules for natural datasets
#


# Only enable downloading of natural sequence and titer data if the user has access to the lab's internal database.
# Otherwise, direct the user how to obtain the correct sequences from GISAID.
if "RETHINK_HOST" in os.environ and "RETHINK_AUTH_KEY" in os.environ:
    localrules: download_sequences, download_all_titers_by_assay, filter_metadata, filter

    rule download_sequences:
        output:
            sequences = "data/natural/{sample}/original_sequences.fasta"
        params:
            fasta_fields = _get_fauna_fields,
            lineage = _get_lineage,
            segment = _get_segment
        conda: "../envs/anaconda.python3.yaml"
        benchmark: "benchmarks/download_sequences_natural_{sample}.txt"
        log: "logs/download_sequences_natural_{sample}.log"
        shell:
            """
            python3 {path_to_fauna}/vdb/download.py \
                --database vdb \
                --virus flu \
                --fasta_fields {params.fasta_fields} \
                --resolve_method split_passage \
                --select locus:{params.segment} lineage:seasonal_{params.lineage} \
                --path data/natural/{wildcards.sample} \
                --fstem original_sequences
            """


    rule download_all_titers_by_assay:
        output:
            titers = DATA_NATURAL_ROOT_PATH + "complete_titers.json"
        conda: "../envs/anaconda.python3.yaml"
        benchmark: "benchmarks/download_all_titers_natural_{sample}.txt"
        log: "logs/download_all_titers_natural_{sample}.log"
        params:
            databases = _get_titer_databases,
            lineage = _get_lineage,
            assay = _get_titer_assay
        shell:
            """
            python3 {path_to_fauna}/tdb/download.py \
                --database {params.databases} \
                --virus flu \
                --subtype {params.lineage} \
                --select assay_type:{params.assay} \
                --path data/natural/{wildcards.sample} \
                --fstem complete \
                --ftype json
            """


    rule download_all_fra_titers:
        output:
            titers = DATA_NATURAL_ROOT_PATH + "complete_fra_titers.json"
        conda: "../envs/anaconda.python3.yaml"
        benchmark: "benchmarks/download_all_fra_titers_natural_{sample}.txt"
        log: "logs/download_all_fra_titers_natural_{sample}.log"
        params:
            databases = _get_titer_databases,
            lineage = _get_lineage,
            assay = "fra"
        shell:
            """
            python3 {path_to_fauna}/tdb/download.py \
                --database {params.databases} \
                --virus flu \
                --subtype {params.lineage} \
                --select assay_type:{params.assay} \
                --path data/natural/{wildcards.sample} \
                --fstem complete_fra \
                --ftype json
            """


    rule get_titers_by_passage:
        input:
            titers = rules.download_all_titers_by_assay.output.titers
        output:
            titers = DATA_NATURAL_ROOT_PATH + "titers.tsv"
        params:
            passage = _get_titer_passage
        benchmark: "benchmarks/get_titers_natural_{sample}.txt"
        log: "logs/get_titers_natural_{sample}.log"
        conda: "../envs/anaconda.python3.yaml"
        shell:
            """
            python3 scripts/get_titers_by_passage.py \
                --titers {input.titers} \
                --passage-type {params.passage} \
                --output {output.titers}
            """


    rule get_fra_titers_by_passage:
        input:
            titers = rules.download_all_fra_titers.output.titers
        output:
            titers = DATA_NATURAL_ROOT_PATH + "fra_titers.tsv"
        params:
            passage = _get_titer_passage
        benchmark: "benchmarks/get_fra_titers_natural_{sample}.txt"
        log: "logs/get_fra_titers_natural_{sample}.log"
        conda: "../envs/anaconda.python3.yaml"
        shell:
            """
            python3 scripts/get_titers_by_passage.py \
                --titers {input.titers} \
                --passage-type {params.passage} \
                --output {output.titers}
            """


    rule parse:
        input:
            sequences = rules.download_sequences.output.sequences
        output:
            sequences = DATA_NATURAL_ROOT_PATH + "sequences.fasta",
            metadata = DATA_NATURAL_ROOT_PATH + "metadata.tsv"
        params:
            fasta_fields = _get_fasta_fields
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
        input:
            metadata = rules.parse.output.metadata,
            sequences = rules.parse.output.sequences
        output:
            sequences = protected(DATA_NATURAL_ROOT_PATH + "filtered_sequences.fasta")
        params:
            min_length = _get_min_sequence_length,
            exclude = _get_outliers
        conda: "../envs/anaconda.python3.yaml"
        benchmark: "benchmarks/filter_natural_{sample}.txt"
        shell:
            """
            augur filter \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --min-length {params.min_length} \
                --exclude {params.exclude} \
                --exclude-where country=? region=? passage=egg \
                --output {output}
            """


    rule filter_metadata:
        input:
            metadata = rules.parse.output.metadata,
        output:
            metadata = DATA_NATURAL_ROOT_PATH + "filtered_metadata.tsv"
        conda: "../envs/anaconda.python3.yaml"
        shell:
            """
            python3 scripts/filter_strains_with_ambiguous_dates.py \
                --metadata {input.metadata} \
                --output {output.metadata}
            """


    rule select_strains:
        input:
            sequences = rules.filter.output.sequences,
            metadata = rules.filter_metadata.output.metadata,
            titers = rules.get_titers_by_passage.output.titers,
            include = _get_required_strains
        output:
            strains = DATA_NATURAL_ROOT_PATH + "strains.txt"
        params:
            viruses_per_month = _get_viruses_per_month,
            lineage = _get_lineage,
            segment = _get_segment,
            start_date = _get_start_date_for_dataset,
            end_date = _get_end_date_for_dataset,
        conda: "../envs/anaconda.python3.yaml"
        benchmark: "benchmarks/select_strains_natural_{sample}.txt"
        log: "logs/select_strains_natural_{sample}.log"
        shell:
            """
            python3 scripts/select_strains.py \
                --sequences {input.sequences} \
                --metadata {input.metadata} \
                --segments {params.segment} \
                --include {input.include} \
                --lineage {params.lineage} \
                --time-interval {params.start_date} {params.end_date} \
                --viruses_per_month {params.viruses_per_month} \
                --titers {input.titers} \
                --output {output.strains}
            """


    rule extract_strain_metadata:
        input:
            strains = _get_strains_for_natural_data,
            metadata = rules.filter_metadata.output.metadata
        output:
            metadata = protected(DATA_NATURAL_ROOT_PATH + "strains_metadata.tsv")
        conda: "../envs/anaconda.python3.yaml"
        shell:
            """
            python3 scripts/filter_metadata_by_strains.py \
                --metadata {input.metadata} \
                --strains {input.strains} \
                --output {output.metadata}
            """
else:
    logger.warning("\nNote that sequences for natural populations must be downloaded manually from GISAID.")
    logger.warning("See instructions for downloading these sequences in the project's documentation:\n")
    logger.warning("  https://github.com/blab/flu-forecasting#download-sequences-for-natural-populations\n")
