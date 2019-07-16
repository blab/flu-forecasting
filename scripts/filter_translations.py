"""Filter translation sequences by date
"""
import argparse
from augur.reconstruct_sequences import load_alignments
from augur.utils import read_node_data
from Bio import AlignIO, SeqIO, Seq, SeqRecord
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", nargs="+", help="amino acid sequences to filter", required=True)
    parser.add_argument("--branch-lengths", help="JSON with branch lengths and internal node dates estimated by TreeTime", required=True)
    parser.add_argument("--min-date", help="minimum date for sequences to emit", required=True)
    parser.add_argument("--output", nargs="+", help="filtered amino acid sequences (one per input)", required=True)
    args = parser.parse_args()

    # Get min date.
    min_date = pd.to_datetime(args.min_date)

    # Load branch lengths.
    node_data = read_node_data(args.branch_lengths)

    # Write alignments to file.
    for alignment_file, output_file in zip(args.alignment, args.output):
        alignments = AlignIO.read(alignment_file, "fasta")
        new_alignments = []
        for alignment in alignments:
            date = pd.to_datetime(node_data["nodes"][alignment.id]["date"])
            if date >= min_date:
                new_alignments.append(alignment)

        SeqIO.write(new_alignments, output_file, "fasta")
