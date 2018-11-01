import argparse
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("parameters")
    parser.add_argument("output")

    args = parser.parse_args()

    fontsize = 18
    params = {
        'axes.labelsize': fontsize,
        'font.size': fontsize,
        'legend.fontsize': 16,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'text.usetex': False,
        'figure.figsize': [8, 6],
        'savefig.dpi': 300,
        'figure.dpi': 120
    }
    plt.rcParams.update(params)

    df = pd.read_table(args.parameters, keep_default_na=False)
    df["build"] = df.apply(lambda row: "%s, %sv" % (row["year_range"], row["viruses"]), axis=1)

    g = sns.catplot(
        y="predictors",
        x="param",
        hue="predictor",
        row="build",
        data=df,
        ci="sd",
        dodge=0.5,
        join=False,
        kind="point",
        height=10,
        aspect=1.33
    )

    # Check for a single axis and if it doesn't exist, we must have multiple axes.
    try:
        axes = [g.ax]
    except AttributeError:
        axes = g.axes.flatten()

    for ax in axes:
        ax.axvline(0, color="#999999", alpha=0.5)

    g.set_axis_labels("Model parameters", "Predictors")
    plt.tight_layout()
    plt.savefig(args.output)
