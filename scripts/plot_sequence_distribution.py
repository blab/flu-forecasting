"""
Plot a histogram of number of sequences in a given tree by month.
"""
import argparse
from augur.utils import json_to_tree
import json
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="tree JSON for auspice")
    parser.add_argument("output", help="sequences histogram plot")

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

    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)

    tip_dates = [
        node.attr["raw_date"]
        for node in tree.find_clades()
        if node.is_terminal() and not "XX" in node.attr["raw_date"]
    ]
    tip_dates_df = pd.DataFrame({"date": pd.to_datetime(pd.Series(tip_dates, None))})

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax = tip_dates_df.resample("6M", on="date").count().plot(ax=ax, kind="barh", legend=False)

    ax.set_title(args.tree)
    ax.set_xlabel("Number of sequences")
    ax.set_ylabel("Date")

    plt.tight_layout()
    plt.savefig(args.output)
