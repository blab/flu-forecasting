import argparse
from augur.frequencies import TreeKdeFrequencies
import json
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("frequencies", help="frequencies JSON")
    parser.add_argument("output", help="frequency trajectories plot")
    parser.add_argument("--min-frequency", type=float, default=0.1, help="clades must reach at least this frequency at some point to be plotted")

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

    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)
        frequencies = TreeKdeFrequencies.from_json(frequencies_json)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    for clade in frequencies.frequencies:
        if frequencies.frequencies[clade].max() >= args.min_frequency:
            ax.plot(frequencies.pivots, frequencies.frequencies[clade])

    ax.set_title(args.frequencies)
    ax.set_xlabel("Date")
    ax.set_ylabel("Frequency")

    plt.savefig(args.output)
