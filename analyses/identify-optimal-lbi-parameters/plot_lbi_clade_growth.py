"""Plot clade growth by LBI.
"""
import argparse
import json
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot clade growth by LBI.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("model", help="JSON for fitness model")
    parser.add_argument("figure", help="figure of clade growth by LBI")
    parser.add_argument("table", help="table of clade growth by LBI")

    args = parser.parse_args()

    with open(args.model, "r") as fh:
        json_model = json.load(fh)

    df = pd.DataFrame(json_model["data"])
    df["observed_growth"] =  df["observed_freq"] / df["initial_freq"]
    df.to_csv(args.table, sep="\t", header=True, index=False)

    g = sns.lmplot("lbi", "observed_growth", data=df, height=6, aspect=1.33, fit_reg=False)
    g.set_axis_labels("LBI", "Observed clade frequency fold change")

    for ax in g.axes.flatten():
        ax.axhline(y=1.0, alpha=0.5, color="black")
        ax.set_xlim(0, 0.6)
        ax.set_ylim(-0.5, 8)

    plt.savefig(args.figure)
    plt.close()
