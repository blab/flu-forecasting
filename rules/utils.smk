#
# Define helper functions.
#

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

def _get_reference(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["reference"]

def _get_pivot_interval(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["pivot_interval"]

def _get_min_date_for_translation_filter(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    min_date = timepoint - pd.DateOffset(years=config["years_for_titer_alignments"])
    return min_date.strftime("%Y-%m-%d")
