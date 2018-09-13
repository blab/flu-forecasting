import argparse
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy.stats import pearsonr
import seaborn as sns
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dist", "augur")
sys.path.append(augur_path)
from base.fitness_model import get_matthews_correlation_coefficient_for_data_frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="model TSV")
    parser.add_argument("output", help="model fold change scatterplot (PDF, PNG, etc.)")
    parser.add_argument("year_range")
    parser.add_argument("viruses")
    parser.add_argument("predictors")
    parser.add_argument("sample")

    args = parser.parse_args()

    frequency_bins = np.arange(0.1, 1, 0.2)
    df = pd.read_table(args.input)
    df["observed_ratio"] = df["observed_freq"] / df["initial_freq"]
    df["predicted_ratio"] = df["predicted_freq"] / df["initial_freq"]
    df["binned_initial_freq"] = pd.cut(df["initial_freq"], bins=frequency_bins)

    correlation_by_bin = []
    mcc_by_bin = []
    for freq, freq_df in df.groupby("binned_initial_freq"):
        correlation_by_bin.append(pearsonr(freq_df["observed_ratio"], freq_df["predicted_ratio"])[0])

        # Calculate Matthew's correlation coefficient
        mcc_by_bin.append(get_matthews_correlation_coefficient_for_data_frame(freq_df))

    g = sns.FacetGrid(df, col="binned_initial_freq", col_wrap=2, size=4)
    g.map(sns.regplot, "observed_ratio", "predicted_ratio", fit_reg=False)
    g.add_legend()

    for i, ax in enumerate(g.axes.flatten()):
        ax.axhline(y=1, color="#999999", alpha=0.7)
        ax.axvline(x=1, color="#999999", alpha=0.7)
        ax.text(
            0.5,
            0.9,
            "Pearson's $R$ = %.2f" % correlation_by_bin[i],
            transform=ax.transAxes,
            horizontalalignment="left",
            verticalalignment="center"
        )
        ax.text(
            0.5,
            0.8,
            "MCC = %.2f" % mcc_by_bin[i],
            transform=ax.transAxes,
            horizontalalignment="left",
            verticalalignment="center"
        )

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle("year: %s, viruses: %s, predictors: %s, sample: %s" % (
        args.year_range,
        args.viruses,
        args.predictors,
        args.sample
    ))
    plt.savefig(args.output)
