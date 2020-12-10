import argparse
from augur.utils import read_node_data, write_json
import Bio.Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import numpy as np
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate consensus sequence for all non-zero frequency strains and the distance of each strain from the resulting consensus.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="node data JSON containing sequences to find a consensus for")
    parser.add_argument("--frequencies", required=True, help="table of strain frequencies at the current timepoint")
    parser.add_argument("--sequence-attribute", default="aa_sequence", help="attribute in node data JSON containing the sequence data to use")
    parser.add_argument("--frequency-attribute", default="kde_frequency", help="attribute in frequency table representing the frequency data to use")
    parser.add_argument("--output", required=True, help="node data JSON with consensus sequence and distances from the consensus per strain")

    args = parser.parse_args()

    # Load sequence data from a node data JSON file.
    node_sequences = read_node_data(args.sequences)

    # Load frequency data.
    frequencies = pd.read_csv(args.frequencies, sep="\t")

    # Select names of strains with non-zero frequencies.
    strains = set(frequencies.query(f"{args.frequency_attribute} > 0.0")["strain"].values)

    # Select sequences for strains with non-zero frequencies.
    sequences = Bio.Align.MultipleSeqAlignment([
        SeqRecord(
            Seq(
                record[args.sequence_attribute],
            ),
            id=strain
        )
        for strain, record in node_sequences["nodes"].items()
        if strain in strains
    ])

    # Output will store the consensus sequence and the distance of each strain
    # to the consensus.
    output = {
        "nodes": {}
    }

    # Calculate the consensus sequence using a majority-rule approach where we
    # take the most common value in each column.
    consensus = "".join(
        Counter(sequences[:, i]).most_common(1)[0][0]
        for i in range(sequences.get_alignment_length())
    )
    output["consensus"] = consensus

    # Calculate the distance of each strain sequence from the consensus.
    consensus_array = np.frombuffer(consensus.encode(), dtype="S1")
    for sequence in sequences:
        sequence_array = np.frombuffer(str(sequence.seq).encode(), dtype="S1")
        distance = int((consensus_array != sequence_array).sum())
        output["nodes"][sequence.id] = {
            "distance_from_consensus": distance
        }

    # Output the results.
    write_json(output, args.output)
