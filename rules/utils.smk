#
# Define helper functions.
#

def _get_sequences_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["sequences"]

def _get_metadata_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["metadata"]

def _get_strains_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["strains"]

def _get_titers_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["titers"]

def _get_start_date_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["start_date"]

def _get_end_date_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["end_date"]

def _get_min_date_for_augur_frequencies_by_wildcards(wildcards):
    return timestamp_to_float(pd.to_datetime(_get_start_date_by_wildcards(wildcards)))

def _get_max_date_for_augur_frequencies_by_wildcards(wildcards):
    return timestamp_to_float(pd.to_datetime(wildcards.timepoint))

def _get_viruses_per_month(wildcards):
    return config["datasets"][wildcards.sample]["viruses_per_month"]

def _get_simulation_seed(wildcards):
    return config["datasets"][wildcards.sample]["seed"]

def _get_fauna_fields(wildcards):
    return config["datasets"][wildcards.sample]["fauna_fields"]

def _get_fasta_fields(wildcards):
    return config["datasets"][wildcards.sample]["fasta_fields"]

def _get_lineage(wildcards):
    return config["datasets"][wildcards.sample]["lineage"]

def _get_segment(wildcards):
    return config["datasets"][wildcards.sample]["segment"]

def _get_titer_databases(wildcards):
    return config["datasets"][wildcards.sample]["titer_databases"]

def _get_titer_assay(wildcards):
    return config["datasets"][wildcards.sample]["titer_assay"]

def _get_titer_passage(wildcards):
    return config["datasets"][wildcards.sample]["titer_passage"]

def _get_min_sequence_length(wildcards):
    return config["datasets"][wildcards.sample]["min_sequence_length"]

def _get_outliers(wildcards):
    return config["datasets"][wildcards.sample]["outliers"]

def _get_required_strains(wildcards):
    return config["datasets"][wildcards.sample]["required_strains"]

def _get_start_date_for_dataset(wildcards):
    return config["datasets"][wildcards.sample]["start_date"]

def _get_end_date_for_dataset(wildcards):
    return config["datasets"][wildcards.sample]["end_date"]

def _get_pivot_interval(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["pivot_interval"]

def _get_min_date_for_translation_filter(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    min_date = timepoint - pd.DateOffset(years=config["years_for_titer_alignments"])
    return min_date.strftime("%Y-%m-%d")
