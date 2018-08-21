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
import datetime
import json
import numpy as np
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
print(augur_path)
sys.path.append(augur_path)
from base.frequencies import KdeFrequencies
from base.io_util import json_to_tree
from base.process import process


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("tree", help="auspice JSON tree")
    parser.add_argument("frequencies", help="JSON with frequencies estimated from the given tree and used to estimate the given parameters")
    parser.add_argument("--narrow-bandwidth", type=float, default=1 / 12.0, help="the bandwidth for the narrow KDE")
    parser.add_argument("--wide-bandwidth", type=float, default=3 / 12.0, help="the bandwidth for the wide KDE")
    parser.add_argument("--proportion-wide", type=float, default=0.2, help="the proportion of the wide bandwidth to use in the KDE mixture model")
    parser.add_argument("--pivot-frequency", type=float, default=1 / 12.0, help="frequency of pivots in fraction of a year")
    parser.add_argument("--start-date", help="the start of the interval to estimate frequencies across")
    parser.add_argument("--end-date", help="the end of the interval to estimate frequencies across")
    parser.add_argument("--include-internal-nodes", action="store_true", help="calculate frequencies for internal nodes as well as tips")

    parser.add_argument("--precision", type=int, default=6, help="number of decimal places to retain in frequency estimates")

    args = parser.parse_args()

    # Load JSON tree.
    with open(args.tree, "r") as json_fh:
         json_tree = json.load(json_fh)

    # Convert JSON tree layout to a Biopython Clade instance.
    tree = json_to_tree(json_tree)

    # Convert start and end dates to floats from time interval format.
    if args.start_date is not None and args.end_date is not None:
        # Convert the string time interval to a datetime instance and then to floats.
        time_interval = [
            datetime.datetime.strptime(time, "%Y-%m-%d")
            for time in (args.end_date, args.start_date)
        ]

        start_date, end_date = process.get_time_interval_as_floats(time_interval)
    else:
        start_date = end_date = None

    # Estimate frequencies.
    frequencies = KdeFrequencies(
        sigma_narrow=args.narrow_bandwidth,
        sigma_wide=args.wide_bandwidth,
        proportion_wide=args.proportion_wide,
        pivot_frequency=args.pivot_frequency,
        start_date=start_date,
        end_date=end_date,
        include_internal_nodes=args.include_internal_nodes
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
        json.dump(json_frequencies, oh, indent=2)
