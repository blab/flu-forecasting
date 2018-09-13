"""Plot model accuracies by LBI parameters.
"""
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Configure matplotlib theme.
fontsize = 14
matplotlib_params = {
    'axes.labelsize': fontsize,
    'font.size': fontsize,
    'legend.fontsize': 12,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'text.usetex': False,
    'figure.figsize': [8, 6],
    'savefig.dpi': 300,
    'figure.dpi': 300,
    'text.usetex': False
}

plt.rcParams.update(matplotlib_params)

# Turn off spines for all plots.
plt.rc("axes.spines", top=False, right=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("results", help="tab-delimited model results for all LBI parameters")
    parser.add_argument("correlation", help="plot of correlations by LBI parameters")
    parser.add_argument("mcc", help="plot of MCCs by LBI parameters")

    args = parser.parse_args()

    results_df = pd.read_table(args.results)

    sns.lmplot(x="tau", y="correlation", hue="time_window", fit_reg=False, data=results_df, aspect=1.5)
    plt.savefig(args.correlation)

    sns.lmplot(x="tau", y="mcc", hue="time_window", fit_reg=False, data=results_df, aspect=1.5)
    plt.savefig(args.mcc)

    print(results_df.to_string())
