#
# Define helper functions.
#
import datetime


def float_to_datestring(time):
    """Convert a floating point date from TreeTime `numeric_date` to a date string
    """
    # Extract the year and remainder from the floating point date.
    year = int(time)
    remainder = time - year

    # Calculate the day of the year (out of 365 + 0.25 for leap years).
    tm_yday = int(remainder * 365.25)
    if tm_yday == 0:
        tm_yday = 1

    # Construct a date object from the year and day of the year.
    date = datetime.datetime.strptime("%s-%s" % (year, tm_yday), "%Y-%j")

    # Build the date string with zero-padded months and days.
    date_string = "%s-%.2i-%.2i" % (date.year, date.month, date.day)

    return date_string

def _get_sample_by_wildcards(wildcards):
    """Look for an explicit dataset sample name for the build sample and default to
    the build sample name when no explicit name is provided.
    """
    try:
        sample = config["builds"][wildcards.type][wildcards.sample].get("dataset_sample", wildcards.sample)
    except AttributeError:
        sample = wildcards.sample

    return sample

def _get_sequences_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["sequences"]

def _get_metadata_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["metadata"]

def _get_complete_metadata_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["complete_metadata"]

def _get_strains_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["strains"]

def _get_titers_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["titers"]

def _get_fra_titers_by_wildcards(wildcards):
    try:
        return config["builds"][wildcards.type][wildcards.sample]["fra_titers"]
    except KeyError:
        print("ERROR: Could not find FRA titers for the %s sample '%s'. Add a 'fra_titers' entry to the config file." % (wildcards.type, wildcards.sample), file=sys.stderr)
        raise

def _get_start_date_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["start_date"]

def _get_end_date_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["end_date"]

def _get_min_date_for_augur_frequencies_by_wildcards(wildcards):
    # Calculate min date for augur frequencies based on the current timepoint
    # minus the maximum number of years back for building trees.
    years_back = _get_years_back_to_build_trees(wildcards)
    offset = pd.DateOffset(years=years_back)
    start_date = pd.to_datetime(_get_start_date_by_wildcards(wildcards))
    timepoint_date = pd.to_datetime(wildcards.timepoint)
    min_date = max(start_date, timepoint_date - offset)

    return timestamp_to_float(min_date)

def _get_min_date_for_diffusion_frequencies_by_wildcards(wildcards):
    # Calculate min date for diffusion frequencies based on the current
    # timepoint minus the maximum number of years back allowed.
    years_back = config["frequencies"]["max_years_for_diffusion"]
    offset = pd.DateOffset(years=years_back)

    start_date = pd.to_datetime(_get_start_date_by_wildcards(wildcards))
    timepoint_date = pd.to_datetime(wildcards.timepoint)
    min_diffusion_date = max(start_date, timepoint_date - offset)

    return timestamp_to_float(min_diffusion_date)

def _get_max_date_for_augur_frequencies_by_wildcards(wildcards):
    return timestamp_to_float(pd.to_datetime(wildcards.timepoint))

def _get_viruses_per_month(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["viruses_per_month"]

def _get_simulation_seed(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["seed"]

def _get_fauna_fields(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["fauna_fields"]

def _get_fasta_fields(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["fasta_fields"]

def _get_lineage(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["lineage"]

def _get_segment(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["segment"]

def _get_titer_databases(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["titer_databases"]

def _get_titer_assay(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["titer_assay"]

def _get_titer_passage(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["titer_passage"]

def _get_min_sequence_length(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["min_sequence_length"]

def _get_outliers(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["outliers"]

def _get_required_strains(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["required_strains"]

def _get_required_strains_argument(wildcards):
    if wildcards.type == "natural":
        required_strains = _get_required_strains(wildcards)
        return "--reference-strains %s" % required_strains
    else:
        # No required strains argument needed for simulated strains.
        return ""

def _get_strains_for_natural_data(wildcards):
    # If the current build has a predefined path to a frozen list of strains to
    # use, return that path. Otherwise, we need to select a new set of strains
    # for a new analysis.
    sample = _get_sample_by_wildcards(wildcards)
    if "frozen_strains" in config["datasets"][sample]:
        return config["datasets"][sample]["frozen_strains"]
    else:
        return rules.select_strains.output.strains

def _get_start_date_for_dataset(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["start_date"]

def _get_end_date_for_dataset(wildcards):
    sample = _get_sample_by_wildcards(wildcards)
    return config["datasets"][sample]["end_date"]

def _get_reference(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["reference"]

def _get_pivot_interval(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["pivot_interval"]

def _get_min_date_for_translation_filter(wildcards):
    timepoint = pd.to_datetime(wildcards.timepoint)
    min_date = timepoint - pd.DateOffset(years=config["years_for_titer_alignments"])
    return min_date.strftime("%Y-%m-%d")

def _get_tree_plots_by_wildcards(wildcards):
    build = config["builds"][wildcards.type][wildcards.sample]
    timepoints = _get_timepoints_for_build_interval(
        build["start_date"],
        build["end_date"],
        build["pivot_interval"],
        build["min_years_per_build"]
    )
    return expand(
        BUILD_PATH.replace("{", "{{").replace("}", "}}") + "timepoints/{timepoint}/tree.pdf",
        timepoint=timepoints
    )

def _get_tip_attributes_by_wildcards(wildcards):
    build = config["builds"][wildcards.type][wildcards.sample]
    timepoints = _get_timepoints_for_build_interval(
        build["start_date"],
        build["end_date"],
        build["pivot_interval"],
        build["min_years_per_build"]
    )
    return expand(
        BUILD_PATH.replace("{", "{{").replace("}", "}}") + "timepoints/{timepoint}/tip_attributes.tsv",
        timepoint=timepoints
    )

def _get_tip_clades_by_wildcards(wildcards):
    build = config["builds"][wildcards.type][wildcards.sample]
    timepoints = _get_timepoints_for_build_interval(
        build["start_date"],
        build["end_date"],
        build["pivot_interval"],
        build["min_years_per_build"]
    )
    return expand(
        BUILD_PATH.replace("{", "{{").replace("}", "}}") + "timepoints/{timepoint}/tips_to_clades.tsv",
        timepoint=timepoints
    )

def _get_final_tree_for_wildcards(wildcards):
    end_date = _get_end_date_by_wildcards(wildcards)
    return BUILD_PATH.format(**wildcards) + "timepoints/%s/tree.nwk" % end_date

def _get_years_back_to_build_trees(wildcards):
    return config["builds"][wildcards.type][wildcards.sample]["years_back_to_build_trees"]

def _get_fitness_model_training_window(wildcards):
    # Prefer a build-specific training window for fitnesses models and default
    # to the global setting if a build-specific value isn't provided.
    return config["builds"][wildcards.type][wildcards.sample].get(
        "training_window",
        config["fitness_model"]["training_window"]
    )

def _get_model_from_validation(wildcards):
    # Convert the path for the best model's minimal JSON to the corresponding
    # path for the model's validation errors.
    path = Path(config["builds"][wildcards.type][wildcards.sample]["best_predictor"])
    return str(path.parent.parent / Path("minimal_models_by_distances") / Path(wildcards.predictors + ".json"))

def _get_model_validation_errors(wildcards):
    # Convert the path for the best model's minimal JSON to the corresponding
    # path for the model's validation errors.
    path = Path(config["builds"][wildcards.type][wildcards.sample]["best_predictor"])
    return str(path.parent.parent / Path("models_by_distances_errors") / Path(wildcards.predictors + ".tsv"))

def _get_delta_months_to_forecast(wildcards):
    return " ".join([str(month) for month in config["fitness_model"]["delta_months"]])

def _get_model_to_test_by_wildcards(wildcards):
    build = config["builds"][wildcards.type][wildcards.sample]
    model_path = rules.extract_minimal_models_by_distances.output.model.format(
        type=wildcards.type,
        sample=build["validation_build"],
        predictors=wildcards.predictors
    )
    return model_path

def _get_validation_sample_by_wildcards(wildcards):
    return config["builds"][wildcards.type][wildcards.sample].get("validation_build", wildcards.sample)

def _get_full_tree_sample_by_wildcards(wildcards):
    if "full_tree_build" not in config["builds"][wildcards.type][wildcards.sample]:
        logger.error("The sample '%s' needs a 'full_tree_build' entry in its build configuration." % wildcards.sample)

    return config["builds"][wildcards.type][wildcards.sample]["full_tree_build"]
