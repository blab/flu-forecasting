import argparse
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


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

    g = sns.FacetGrid(df, col="binned_initial_freq", col_wrap=2, size=4)
    g.map(sns.regplot, "observed_ratio", "predicted_ratio", fit_reg=False)
    g.add_legend()

    for ax in g.axes.flatten():
        ax.axhline(y=1, color="#999999", alpha=0.7)
        ax.axvline(x=1, color="#999999", alpha=0.7)

    plt.subplots_adjust(top=0.9)
    g.fig.suptitle("year: %s, viruses: %s, predictors: %s, sample: %s" % (
        args.year_range,
        args.viruses,
        args.predictors,
        args.sample
    ))
    plt.savefig(args.output)
