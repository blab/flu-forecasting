"""Partition a given list of strains by their metadata into the given timepoints.
"""
import argparse
from augur.utils import get_numerical_dates, read_metadata
import numpy as np
import pandas as pd
from treetime.utils import numeric_date

from select_strains import read_strain_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Partition strains into timepoints",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("metadata", help="tab-delimited metadata with columns for strain and date")
    parser.add_argument("timepoint", help="date for which strains should be partitioned")
    parser.add_argument("output", help="text file into which strains should be written for the given timepoint")
    parser.add_argument("--years-back", type=int, help="Number of years prior to the given timepoint to limit strains to")
    parser.add_argument("--additional-years-back-for-references", type=int, default=5, help="Additional number of years prior to the given timepoint to allow reference strains")
    parser.add_argument("--reference-strains", help="text file containing list of reference strains that should be included from the original strains even if they were sampled prior to the minimum date determined by the requested number of years before the given timepoint")
    args = parser.parse_args()

    # Convert date string to a datetime instance.
    timepoint = pd.to_datetime(args.timepoint)
    numeric_timepoint = np.around(numeric_date(timepoint), 2)

    # Load metadata with strain names and dates.
    metadata, columns = read_metadata(args.metadata)

    # Convert string dates with potential ambiguity (e.g., 2010-05-XX) into
    # floating point dates.
    dates = get_numerical_dates(metadata, fmt="%Y-%m-%d")

    # Setup reference strains.
    if args.reference_strains:
        reference_strains = read_strain_list(args.reference_strains)
    else:
        reference_strains = []

    # If a given number of years back has been requested, determine what the
    # earliest date to accept for strains is.
    if args.years_back is not None:
        earliest_timepoint = timepoint - pd.DateOffset(years=args.years_back)
        numeric_earliest_timepoint = np.around(numeric_date(earliest_timepoint), 2)

        # If reference strains are provided, calculate the earliest date to
        # accept those strains.
        if len(reference_strains) > 0:
            earliest_reference_timepoint = earliest_timepoint - pd.DateOffset(years=args.additional_years_back_for_references)
            numeric_earliest_reference_timepoint = np.around(numeric_date(earliest_reference_timepoint), 2)

    # Find strains sampled prior to the current timepoint. Strains may have
    # multiple numerical dates, so we filter on the latest (maximum) observed
    # date per strain. If a requested number of years back is provided, use the
    # corresponding earliest dates for non-reference and reference strains to
    # determine whether they are included in the current timepoint.
    timepoint_strains = []
    for strain, strain_dates in dates.items():
        strain_date = np.max(strain_dates)
        if (strain_date <= numeric_timepoint and
            ((args.years_back is None) or
             (strain_date >= numeric_earliest_timepoint) or
             (strain in reference_strains and strain_date >= numeric_earliest_reference_timepoint))):
            timepoint_strains.append(strain)

    timepoint_strains = sorted(timepoint_strains)

    # Write sorted list of strains to disk.
    with open(args.output, "w") as oh:
        for strain in timepoint_strains:
            oh.write(f"{strain}\n")
