"""Select sequences for strains that have titer measurements by a given timepoint.
"""
import argparse
from augur.titer_model import TiterCollection
from augur.utils import read_metadata
import Bio.SeqIO
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Select sequences for strains that have titer measurements by a given timepoint",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--titers", required=True, help="table of titer measurements without a header")
    parser.add_argument("--metadata", required=True, help="strain metadata with a date field")
    parser.add_argument("--sequences", required=True, help="FASTA file of strain sequences")
    parser.add_argument("--timepoint", required=True, help="latest date in YYYY-MM-DD format to limit selected sequences to")
    parser.add_argument("--years-back", required=True, type=int, help="number of years prior to the given timepoint to include strains used as test viruses in titer measurements")
    parser.add_argument('--output', help="name of the file to write selected strains to")

    args = parser.parse_args()

    # Setup dates.
    latest_timepoint = pd.to_datetime(args.timepoint)
    offset = pd.DateOffset(years=args.years_back)
    earliest_timepoint = latest_timepoint - offset

    # Load titers.
    measurements, strains, sources = TiterCollection.load_from_file(args.titers)
    test_strains = set([key[0] for key in measurements.keys()])
    reference_strains = set([key[1][0] for key in measurements.keys()])

    # Read in metadata, parse numeric dates
    metadata = pd.read_csv(args.metadata, sep="\t")
    print("Loaded metadata for %i strains" % metadata.shape[0])

    # Drop strains with ambiguous dates and egg passaged status.
    metadata = metadata[
        (~metadata["date"].str.contains("XX")) &
        (~metadata["strain"].str.contains("-egg"))
    ].copy()
    metadata["date"] = pd.to_datetime(metadata["date"])
    print("Retained metadata for %i strains with unambiguous dates and no egg passages" % metadata.shape[0])

    # Annotate strains by titer test or reference status.
    metadata["is_test_virus"] = metadata["strain"].isin(test_strains)
    metadata["is_reference_virus"] = metadata["strain"].isin(reference_strains)

    # Test viruses can only be included if they were sampled prior to the
    # timepoint and after the given number of years back from the
    # timepoint. This limits recurrent mutations that may have different
    # antigenic properties based on their historical context.
    test_strain_indices = (
        (metadata["date"] >= earliest_timepoint) &
        (metadata["date"] <= latest_timepoint) &
        (metadata["is_test_virus"])
    )
    print("Found %i matching test strains" % test_strain_indices.sum())

    # References can be included any time before the timepoint. These sequences
    # are required to allow the model to have access to all available titer
    # measurements.
    reference_strain_indices = (
        (metadata["date"] <= latest_timepoint) &
        (metadata["is_reference_virus"])
    )
    print("Found %i matching reference strains" % reference_strain_indices.sum())

    # Filter strains by date and titer status.
    filtered_metadata = metadata[
        (test_strain_indices) | (reference_strain_indices)
    ].copy()

    print("Found %i strains with metadata" % filtered_metadata.shape[0])

    # Extract sequences for the selected strains.
    selected_strains = filtered_metadata["strain"].values
    sequences = Bio.SeqIO.parse(args.sequences, "fasta")
    records = [sequence for sequence in sequences if sequence.id in selected_strains]
    print("Found %i sequences matching selected strains" % len(records))

    Bio.SeqIO.write(records, args.output, "fasta")
