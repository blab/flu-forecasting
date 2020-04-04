import argparse
import Bio.SeqIO
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", help="simulated sequences that have already been filtered")
    parser.add_argument("--metadata", help="original metadata table for simulated sequences")
    parser.add_argument("--output", help="filtered metadata where only samples present in the given sequences are included")

    args = parser.parse_args()

    # Get a list of all samples that passed the sequence filtering step.
    sequences = Bio.SeqIO.parse(args.sequences, "fasta")
    sample_ids = [sequence.id for sequence in sequences]

    # Load all metadata.
    metadata = pd.read_csv(args.metadata, sep="\t")
    filtered_metadata = metadata[metadata["strain"].isin(sample_ids)].copy()

    # Save only the metadata records that have entries in the filtered sequences.
    filtered_metadata.to_csv(args.output, sep="\t", header=True, index=False)
