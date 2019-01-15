"""Partition a given list of strains by their metadata into the given timepoints.
"""
import argparse
from augur.utils import get_numerical_dates, read_metadata
import numpy as np
import pandas as pd
from treetime.utils import numeric_date


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Partition strains into timepoints",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("metadata", help="tab-delimited metadata with columns for strain and date")
    parser.add_argument("timepoint", help="date for which strains should be partitioned")
    parser.add_argument("output", help="text file into which strains should be written for the given timepoint")
    args = parser.parse_args()

    # Convert date string to a datetime instance.
    timepoint = pd.to_datetime(args.timepoint)
    numeric_timepoint = numeric_date(timepoint)

    # Load metadata with strain names and dates.
    metadata, columns = read_metadata(args.metadata)

    # Convert string dates with potential ambiguity (e.g., 2010-05-XX) into
    # floating point dates.
    dates = get_numerical_dates(metadata, fmt="%Y-%m-%d")

    # Find strains sampled prior to the current timepoint. Strains may have
    # multiple numerical dates, so we filter on the latest (maximum) observed
    # date per strain.
    timepoint_strains = sorted([
        strain
        for strain, strain_dates in dates.items()
        if np.max(strain_dates) <= numeric_timepoint
    ])

    # Write sorted list of strains to disk.
    with open(args.output, "w") as oh:
        for strain in timepoint_strains:
            oh.write(f"{strain}\n")
