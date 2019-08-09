"""Normalize fitness by timepoint frequencies for samples in simulated populations.
"""
import argparse
from augur.utils import write_json
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Normalize fitness by timepoint frequencies for samples in simulated populations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="file with metadata associated with viral sequences, one for each segment")
    parser.add_argument("--frequencies-table", required=True, help="frequencies table for the current timepoint")
    parser.add_argument("--frequency-method", required=True, choices=["kde", "diffusion"], help="method used to estimate frequencies")
    parser.add_argument("--output", required=True, help="JSON of normalized fitness per sample")

    args = parser.parse_args()

    # Load metadata.
    metadata = pd.read_csv(args.metadata, sep="\t")

    # Load frequencies.
    frequencies = pd.read_csv(args.frequencies_table, sep="\t")

    # Filter samples to those with nonzero frequencies at the current timepoint.
    nonzero_frequencies = frequencies[frequencies["%s_frequency" % args.frequency_method] > 0].copy()

    # Merge extent sample frequencies with metadata containing fitnesses.
    nonzero_metadata = nonzero_frequencies.merge(
        metadata,
        on="strain"
    )

    # Normalize fitness by maximum fitness.
    nonzero_metadata["normalized_fitness"] = nonzero_metadata["fitness"] / nonzero_metadata["fitness"].max()

    # Prepare dictionary of normalized fitnesses by sample.
    normalized_fitness = {
        strain: {"normalized_fitness": fitness}
        for strain, fitness in nonzero_metadata.loc[:, ["strain", "normalized_fitness"]].values
    }

    print("Raw fitness: %.2f +/- %.2f" % (nonzero_metadata["fitness"].mean(),
                                          nonzero_metadata["fitness"].std()))
    print("Normalized fitness: %.2f +/- %.2f" % (nonzero_metadata["normalized_fitness"].mean(),
                                                 nonzero_metadata["normalized_fitness"].std()))

    # Save normalized fitness as a node data JSON.
    write_json({"nodes": normalized_fitness}, args.output)
