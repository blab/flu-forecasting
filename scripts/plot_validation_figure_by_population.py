"""
Produce all validation/test figures for all populations in a single notebook.
"""
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
import statsmodels.api as sm


np.random.seed(314159)

PLOT_THEME_ATTRIBUTES = {
    "axes.labelsize": 14,
    "font.size": 18,
    "legend.fontsize": 12,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "figure.figsize": [6.0, 4.0],
    "savefig.dpi": 200,
    "figure.dpi": 200,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.usetex": False
}


def matthews_correlation_coefficient(tp, tn, fp, fn):
    """Return Matthews correlation coefficient for values from a confusion matrix.
    Implementation is based on the definition from wikipedia:

    https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
    """
    numerator = (tp * tn) - (fp * fn)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
            denominator = 1

    return float(numerator) / denominator


def get_matthews_correlation_coefficient_for_data_frame(freq_df, return_confusion_matrix=False):
        """Calculate Matthew's correlation coefficient from a given pandas data frame
        with columns for initial, observed, and predicted frequencies.
        """
        observed_growth = (freq_df["frequency_final"] > freq_df["frequency"])
        predicted_growth = (freq_df["projected_frequency"] > freq_df["frequency"])
        true_positives = ((observed_growth) & (predicted_growth)).sum()
        false_positives= ((~observed_growth) & (predicted_growth)).sum()

        observed_decline = (freq_df["frequency_final"] < freq_df["frequency"])
        predicted_decline = (freq_df["projected_frequency"] < freq_df["frequency"])
        true_negatives = ((observed_decline) & (predicted_decline)).sum()
        false_negatives = ((~observed_decline) & (predicted_decline)).sum()

        mcc = matthews_correlation_coefficient(
            true_positives,
            true_negatives,
            false_positives,
            false_negatives
        )

        if return_confusion_matrix:
            confusion_matrix = {
                "tp": true_positives,
                "tn": true_negatives,
                "fp": false_positives,
                "fn": false_negatives
            }

            return mcc, confusion_matrix
        else:
            return mcc


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints with standardized predictors and weighted distances to the future")
    parser.add_argument("--tips-to-clades", required=True, help="tab-delimited file of all clades per tip and timepoint from a single tree that includes all tips in the given tip attributes table")
    parser.add_argument("--forecasts", required=True, help="table of forecasts for the given tips")
    parser.add_argument("--model-errors", required=True, help="annotated validation errors for the model used to make the given forecasts")
    parser.add_argument("--bootstrap-samples", type=int, default=100, help="number of bootstrap samples to generate for confidence intervals around absolute forecast errors")
    parser.add_argument("--population", help="the population being analyzed (e.g., simulated or natural)")
    parser.add_argument("--sample", help="sample name for population being analyzed")
    parser.add_argument("--predictors", help="predictors being analyzed")
    parser.add_argument("--output", required=True, help="validation figure")
    parser.add_argument("--output-clades-table", help="table of clade frequencies used in left panels")
    parser.add_argument("--output-ranks-table", help="table of strain ranks used in right panels")

    args = parser.parse_args()

    # Define constants for frequency analyses below.
    min_clade_frequency = 0.15
    precision = 4
    pseudofrequency = 0.001
    number_of_bootstrap_samples = args.bootstrap_samples

    sns.set_style("white")
    mpl.rcParams.update(PLOT_THEME_ATTRIBUTES)

    # Load validation errors for the model used to produce the given forecasts
    # table. These errors are used to identify the first validation timepoint.
    model_errors = pd.read_csv(
        args.model_errors,
        sep="\t",
        parse_dates=["validation_timepoint"]
    )
    first_validation_timepoint = model_errors["validation_timepoint"].min().strftime("%Y-%m-%d")

    # Load tip attributes to be associated with clades and used to calculate
    # clade frequencies.
    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        parse_dates=["timepoint"],
        usecols=["strain", "timepoint", "frequency", "aa_sequence"]
    )
    tips = tips.query("timepoint >= '%s'" % first_validation_timepoint).copy()
    distinct_tips_with_sequence = tips.groupby(["timepoint", "aa_sequence"]).first().reset_index()

    # Load mapping of tips to clades based on a single tree that included all of
    # the tips in the given tip attributes table.
    tips_to_clades = pd.read_csv(
        args.tips_to_clades,
        sep="\t",
        usecols=["tip", "clade_membership", "depth"]
    )
    tips_to_clades = tips_to_clades.rename(columns={"tip": "strain"})

    # Load forecasts for all tips by the model associated with the given model
    # errors. First, load only a subset of the forecast information to simplify
    # downstream data frames.
    forecasts = pd.read_csv(
        args.forecasts,
        sep="\t",
        parse_dates=["timepoint"],
        usecols=["timepoint", "strain", "frequency", "projected_frequency"]
    )

    # Next, load the complete forecasts data frame for ranking of estimated and
    # observed closest strains.
    full_forecasts = pd.read_csv(
        args.forecasts,
        sep="\t",
        parse_dates=["timepoint", "future_timepoint"]
    )
    full_forecasts = full_forecasts.query("timepoint >= '%s'" % first_validation_timepoint).copy()

    # Map tip attributes to all corresponding clades.
    clade_tip_initial_frequencies = tips_to_clades.merge(
        tips,
        on=["strain"]
    )
    clade_tip_initial_frequencies["future_timepoint"] = clade_tip_initial_frequencies["timepoint"] + pd.DateOffset(months=12)

    # Calculate the initial frequency of each clade per timepoint.
    initial_clade_frequencies = clade_tip_initial_frequencies.groupby([
        "timepoint", "future_timepoint", "clade_membership"
    ])["frequency"].sum().reset_index()

    # Merge clade frequencies between adjacent years.
    initial_and_observed_clade_frequencies = initial_clade_frequencies.merge(
        initial_clade_frequencies,
        left_on=["future_timepoint", "clade_membership"],
        right_on=["timepoint", "clade_membership"],
        suffixes=["", "_final"]
    ).groupby(["timepoint", "clade_membership", "frequency"])["frequency_final"].sum().reset_index()

    # Select clades with an initial frequency above the defined threshold.
    large_clades = initial_and_observed_clade_frequencies.query("frequency > %s" % min_clade_frequency).copy()

    # Find estimated future frequencies of large clades.
    clade_tip_estimated_frequencies = tips_to_clades.merge(
        forecasts,
        on=["strain"]
    )
    estimated_clade_frequencies = clade_tip_estimated_frequencies.groupby(
        ["timepoint", "clade_membership"]
    ).aggregate({"projected_frequency": "sum"}).reset_index()

    # Annotate initial and observed clade frequencies with the estimated future
    # values.
    complete_clade_frequencies = large_clades.merge(
        estimated_clade_frequencies,
        on=["timepoint", "clade_membership"],
        suffixes=["", "_other"]
    )

    # Reduce precision of frequency estimates to a reasonable value and
    # eliminate entries where the clade frequency did not change between the
    # initial and final timepoints (these are primarily clades that have already
    # fixed at 100%).
    complete_clade_frequencies = np.round(complete_clade_frequencies, 2)
    complete_clade_frequencies = complete_clade_frequencies.query("frequency != frequency_final").copy()

    # Calculate accuracy of growth and decline classifications.
    mcc, confusion_matrix = get_matthews_correlation_coefficient_for_data_frame(complete_clade_frequencies, True)
    growth_accuracy = confusion_matrix["tp"] / float(confusion_matrix["tp"] + confusion_matrix["fp"])
    decline_accuracy = confusion_matrix["tn"] / float(confusion_matrix["tn"] + confusion_matrix["fn"])

    # Calculate the observed and estimated log growth rates for all clades.
    complete_clade_frequencies["log_observed_growth_rate"] = (
        np.log10((complete_clade_frequencies["frequency_final"] + pseudofrequency) / (complete_clade_frequencies["frequency"] + pseudofrequency))
    )
    complete_clade_frequencies["log_estimated_growth_rate"] = (
        np.log10((complete_clade_frequencies["projected_frequency"] + pseudofrequency) / (complete_clade_frequencies["frequency"] + pseudofrequency))
    )

    # Calculate the bounds for the clade growth rate display based on values in
    # observed and estimated rates.
    log_lower_limit = complete_clade_frequencies.loc[:, ["log_observed_growth_rate", "log_estimated_growth_rate"]].min().min() - 0.1
    log_upper_limit = np.ceil(complete_clade_frequencies.loc[:, ["log_observed_growth_rate", "log_estimated_growth_rate"]].max().max()) + 0.1

    # Calculate the Pearson's correlation between observed and estimated log
    # growth rates.
    r, p = pearsonr(
        complete_clade_frequencies["log_observed_growth_rate"],
        complete_clade_frequencies["log_estimated_growth_rate"]
    )

    # Use observed forecasting errors to inspect the accuracy of one-year
    # lookaheads based on the initial frequency of each clade.
    complete_clade_frequencies["clade_error"] = complete_clade_frequencies["frequency_final"] - complete_clade_frequencies["projected_frequency"]
    complete_clade_frequencies["absolute_clade_error"] = np.abs(complete_clade_frequencies["clade_error"])

    # Estimate uncertainty of the mean absolute clade error by initial clade
    # frequency with LOESS fits to bootstraps from the complete data frame.
    bootstrap_samples = []
    for i in range(number_of_bootstrap_samples):
        complete_clade_frequencies_sample = complete_clade_frequencies.sample(frac=1.0, replace=True).copy()
        z = sm.nonparametric.lowess(
            complete_clade_frequencies_sample["absolute_clade_error"].values * 100,
            complete_clade_frequencies_sample["frequency"].values * 100
        )

        # Track both the initial frequency and the LOESS fits for each bootstrap
        # sample. This ensures that the summary statistics calculated downstream
        # per initial frequency are based on the correct LOESS values.
        bootstrap_samples.append(
            pd.DataFrame({
                "initial_frequency": z[:, 0],
                "loess": z[:, 1]}
            )
        )

    bootstrap_df = pd.concat(bootstrap_samples)

    # Calculate the mean and 95% CIs from bootstraps.
    bootstrap_summary = bootstrap_df.groupby("initial_frequency")["loess"].agg(
        lower=lambda group: np.percentile(group, 2.5),
        mean=np.mean,
        upper=lambda group: np.percentile(group, 97.5)
    ).reset_index()

    initial_frequency = bootstrap_summary["initial_frequency"].values
    mean_lowess_fit = bootstrap_summary["mean"].values
    upper_lowess_fit = bootstrap_summary["upper"].values
    lower_lowess_fit = bootstrap_summary["lower"].values

    # For each timepoint, calculate the percentile rank of each strain based on
    # both its observed and estimated distance to the future.
    sorted_df = full_forecasts.dropna().sort_values(
        ["timepoint"]
    ).copy()

    # Filter sorted records by strains with distinct amino acid sequences.
    sorted_df = sorted_df.merge(
        distinct_tips_with_sequence,
        on=["timepoint", "strain"]
    )

    # First, calculate the rank per strain by observed distance to the future.
    sorted_df["timepoint_rank"] = sorted_df.groupby("timepoint")["weighted_distance_to_future"].rank(pct=True)

    # Then, calculate the rank by estimated distance to the future.
    sorted_df["timepoint_estimated_rank"] = sorted_df.groupby("timepoint")["y"].rank(pct=True)

    # Calculate the Spearman correlation of ranks, to get a measure of the model
    # fit.
    rank_rho, rank_p = spearmanr(
        sorted_df["timepoint_rank"],
        sorted_df["timepoint_estimated_rank"]
    )

    # Select the observed rank of the estimated closest strain to the future per
    # timepoint.
    best_fitness_rank_by_timepoint_df = sorted_df.sort_values(
        ["timepoint", "y"],
        ascending=True
    ).groupby("timepoint")["timepoint_rank"].first().reset_index()

    #
    # Summarize model fit by clade frequencies and strain ranks.
    #

    fig = plt.figure(figsize=(10, 10), facecolor='w')
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1], wspace=0.1)

    ticks = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    #
    # Top-left: Clade growth rate correlations
    #

    clade_ax = fig.add_subplot(gs[0])
    clade_ax.plot(
        complete_clade_frequencies["log_observed_growth_rate"],
        complete_clade_frequencies["log_estimated_growth_rate"],
        "o",
        alpha=0.4
    )

    clade_ax.axhline(color="#cccccc", zorder=-5)
    clade_ax.axvline(color="#cccccc", zorder=-5)

    if p < 0.001:
        p_value = "$p value$ < 0.001"
    else:
        p_value = "$p$ = %.3f" % p

    clade_ax.text(
        0.02,
        0.15,
        "Growth accuracy = %.2f\nDecline accuracy = %.2f\nPearson $R^2$ = %.2f\nN = %s" % (
            growth_accuracy,
            decline_accuracy,
            r ** 2,
            complete_clade_frequencies.shape[0]
        ),
        fontsize=12,
        horizontalalignment="left",
        verticalalignment="center",
        transform=clade_ax.transAxes
    )

    clade_ax.set_xlabel("Observed $log_{10}$ fold change")
    clade_ax.set_ylabel("Estimated $log_{10}$ fold change")

    growth_rate_ticks = np.arange(-6, 4, 1)
    clade_ax.set_xticks(growth_rate_ticks)
    clade_ax.set_yticks(growth_rate_ticks)

    clade_ax.set_xlim(log_lower_limit, log_upper_limit)
    clade_ax.set_ylim(log_lower_limit, log_upper_limit)
    clade_ax.set_aspect("equal")

    #
    # Top-right: Estimated closest strain to the future ranking
    #

    rank_ax = fig.add_subplot(gs[1])

    median_best_rank = best_fitness_rank_by_timepoint_df["timepoint_rank"].median()

    rank_ax.hist(best_fitness_rank_by_timepoint_df["timepoint_rank"], bins=np.arange(0, 1.01, 0.05), label=None)
    rank_ax.axvline(
        median_best_rank,
        color="orange",
        label="median = %i%%" % round(median_best_rank * 100, 0)
    )
    rank_ax.set_xticks(ticks)
    rank_ax.set_xticklabels(['{:3.0f}%'.format(x*100) for x in ticks])
    rank_ax.set_xlim(0, 1)

    rank_ax.legend(
        frameon=False
    )
    rank_ax.set_xlabel("Percentile rank by distance\nfor estimated closest strain")
    rank_ax.set_ylabel("Number of timepoints")

    #
    # Bottom-left: Absolute clade forecast errors with uncertainty.
    #

    forecast_error_ax = fig.add_subplot(gs[2])
    forecast_error_ax.plot(
        complete_clade_frequencies["frequency"].values * 100,
        complete_clade_frequencies["absolute_clade_error"].values * 100,
        "o",
        alpha=0.2
    )

    forecast_error_ax.fill_between(
        initial_frequency,
        lower_lowess_fit,
        upper_lowess_fit,
        alpha=0.1,
        color="black"
    )
    forecast_error_ax.plot(
        initial_frequency,
        mean_lowess_fit,
        alpha=0.75,
        color="black"
    )

    forecast_error_ax.set_xlabel("Initial clade frequency")
    forecast_error_ax.set_ylabel("Absolute forecast error")

    forecast_error_ax.set_xticks(ticks * 100)
    forecast_error_ax.set_yticks(ticks * 100)
    forecast_error_ax.set_xticklabels(['{:3.0f}%'.format(x * 100) for x in ticks])
    forecast_error_ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in ticks])

    forecast_error_ax.set_aspect("equal")

    #
    # Bottom-right: Observed vs. estimated percentile rank for all strains at all timepoints.
    #

    all_rank_ax = fig.add_subplot(gs[3])

    if rank_p < 0.001:
        rank_p_value = "$p$ < 0.001"
    else:
        rank_p_value = "$p$ = %.3f" % rank_p

    all_rank_ax.plot(
        sorted_df["timepoint_rank"],
        sorted_df["timepoint_estimated_rank"],
        "o",
        alpha=0.05
    )

    all_rank_ax.text(
        0.45,
        0.05,
        "Spearman $\\rho^2$ = %.2f" % (rank_rho ** 2,),
        fontsize=12,
        horizontalalignment="left",
        verticalalignment="center",
        transform=all_rank_ax.transAxes
    )

    all_rank_ax.set_xticks(ticks)
    all_rank_ax.set_yticks(ticks)
    all_rank_ax.set_xticklabels(['{:3.0f}%'.format(x * 100) for x in ticks])
    all_rank_ax.set_yticklabels(['{:3.0f}%'.format(x * 100) for x in ticks])

    all_rank_ax.set_xlabel("Observed percentile rank")
    all_rank_ax.set_ylabel("Estimated percentile rank")
    all_rank_ax.set_aspect("equal")

    # Annotate panel labels.
    panel_labels_dict = {
        "weight": "bold",
        "size": 14
    }
    plt.figtext(0.0, 0.97, "A", **panel_labels_dict)
    plt.figtext(0.5, 0.97, "B", **panel_labels_dict)
    plt.figtext(0.0, 0.47, "C", **panel_labels_dict)
    plt.figtext(0.5, 0.47, "D", **panel_labels_dict)

    gs.tight_layout(fig)
    plt.savefig(args.output)

    timepoints_better_than_20th_percentile = (best_fitness_rank_by_timepoint_df["timepoint_rank"] <= 0.2).sum()
    total_timepoints = best_fitness_rank_by_timepoint_df.shape[0]
    print(
        "Estimated strain was in the top 20th percentile at %s of %s (%s%%) timepoints" % (
            timepoints_better_than_20th_percentile,
            total_timepoints,
            int(np.round((timepoints_better_than_20th_percentile / float(total_timepoints)) * 100))
        )
    )

    if args.output_clades_table:
        complete_clade_frequencies = complete_clade_frequencies.rename(columns={
            "frequency": "initial_frequency",
            "frequency_final": "observed_future_frequency",
            "projected_frequency": "estimated_future_frequency"
        })
        complete_clade_frequencies["population"] = args.population
        complete_clade_frequencies["predictors"] = args.predictors
        complete_clade_frequencies["error_type"] = "test" if "test" in args.sample else "validation"

        complete_clade_frequencies.to_csv(
            args.output_clades_table,
            sep="\t",
            header=True,
            index=False
        )

    if args.output_ranks_table:
        sorted_df["observed_distance_to_future"] = sorted_df["weighted_distance_to_future"]
        sorted_df["estimated_distance_to_future"] = sorted_df["y"]
        sorted_df["observed_rank"] = sorted_df["timepoint_rank"]
        sorted_df["estimated_rank"] = sorted_df["timepoint_estimated_rank"]

        sorted_df["population"] = args.population
        sorted_df["sample"] = args.sample
        sorted_df["predictors"] = args.predictors
        sorted_df["error_type"] = "test" if "test" in args.sample else "validation"
        sorted_df = np.around(sorted_df, 2)

        sorted_df.to_csv(
            args.output_ranks_table,
            sep="\t",
            header=True,
            index=False,
            columns=[
                "population",
                "error_type",
                "predictors",
                "timepoint",
                "strain",
                "observed_distance_to_future",
                "estimated_distance_to_future",
                "observed_rank",
                "estimated_rank"
            ]
        )
