"""
Estimate frequencies for a given tree and save the output in a JSON file.

Usage:

# Use defaults.
python frequencies.py tree.json frequencies.json

# Specify custom parameters.
python frequencies.py tree.json frequencies.json \
    --narrow-bandwidth 0.25 \
    --start-date 2006-10-01 --end-date 2018-04-01
"""
import argparse
from augur.utils import read_metadata, get_numerical_dates
from augur.frequencies import TreeKdeFrequencies
import Bio
import Bio.Phylo
import datetime
import json
import numpy as np
import os
import sys


def get_time_interval_as_floats(time_interval):
    """
    Converts the given datetime interval to start and end floats.

    Returns:
        start_date (float): the start of the given time interval
        end_date (float): the end of the given time interval
    """
    start_date = time_interval[1].year + (time_interval[1].month - 1) / 12.0
    end_date = time_interval[0].year + (time_interval[0].month - 1) / 12.0
    return start_date, end_date


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("tree", help="Newick tree")
    parser.add_argument("metadata", help="tab-delimited metadata for tips in the given tree including a date field")
    parser.add_argument("frequencies", help="JSON with frequencies estimated from the given tree and used to estimate the given parameters")
    parser.add_argument("--narrow-bandwidth", type=float, default=1 / 12.0, help="the bandwidth for the narrow KDE")
    parser.add_argument("--wide-bandwidth", type=float, default=3 / 12.0, help="the bandwidth for the wide KDE")
    parser.add_argument("--proportion-wide", type=float, default=0.2, help="the proportion of the wide bandwidth to use in the KDE mixture model")
    parser.add_argument("--pivot-frequency", type=int, default=1, help="number of months between pivots")
    parser.add_argument("--start-date", help="the start of the interval to estimate frequencies across")
    parser.add_argument("--end-date", help="the end of the interval to estimate frequencies across")
    parser.add_argument("--include-internal-nodes", action="store_true", help="calculate frequencies for internal nodes as well as tips")
    parser.add_argument("--weights", help="a dictionary of key/value mappings in JSON format used to weight tip frequencies")
    parser.add_argument("--weights-attribute", help="name of the attribute on each tip whose values map to the given weights dictionary")

    parser.add_argument("--precision", type=int, default=6, help="number of decimal places to retain in frequency estimates")
    parser.add_argument("--censored", action="store_true", help="calculate censored frequencies at each pivot")

    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load metadata.
    metadata, columns = read_metadata(args.metadata)
    dates = get_numerical_dates(metadata, fmt='%Y-%m-%d')

    # Annotate tree with dates and other metadata.
    for tip in tree.find_clades(terminal=True):
        tip.attr = {"num_date": np.mean(dates[tip.name])}

        # Annotate tips with metadata to enable filtering and weighting of
        # frequencies by metadata attributes.
        for key, value in metadata[tip.name].items():
            tip.attr[key] = value

    # Convert start and end dates to floats from time interval format.
    if args.start_date is not None and args.end_date is not None:
        # Convert the string time interval to a datetime instance and then to floats.
        time_interval = [
            datetime.datetime.strptime(time, "%Y-%m-%d")
            for time in (args.end_date, args.start_date)
        ]

        start_date, end_date = get_time_interval_as_floats(time_interval)
    else:
        start_date = end_date = None

    # Load weights if they have been provided.
    if args.weights:
        with open(args.weights, "r") as fh:
            weights = json.load(fh)

        weights_attribute = args.weights_attribute
    else:
        weights = None
        weights_attribute = None

    # Estimate frequencies.
    frequencies = TreeKdeFrequencies(
        sigma_narrow=args.narrow_bandwidth,
        sigma_wide=args.wide_bandwidth,
        proportion_wide=args.proportion_wide,
        pivot_frequency=args.pivot_frequency,
        start_date=start_date,
        end_date=end_date,
        weights=weights,
        weights_attribute=weights_attribute,
        include_internal_nodes=args.include_internal_nodes,
        censored=args.censored
    )
    frequencies.estimate(tree)

    # Export frequencies to JSON.
    json_frequencies = frequencies.to_json()

    # Set precision of frequency estimates.
    for clade in json_frequencies["data"]["frequencies"]:
        json_frequencies["data"]["frequencies"][clade] = np.around(
            np.array(
                json_frequencies["data"]["frequencies"][clade]
            ),
            args.precision
        ).tolist()

    with open(args.frequencies, "w") as oh:
        json.dump(json_frequencies, oh, indent=1, sort_keys=True)
