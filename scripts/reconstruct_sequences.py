"""Calculate LBI for a given tree and one or more sets of parameters.
"""
import argparse
import Bio
import Bio.Phylo
import json
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.io_util import json_to_tree, reconstruct_sequences_from_mutations


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick file")
    parser.add_argument("nucleotide_mutations", help="JSON of nucleotide sequences and mutations per node")
    parser.add_argument("amino_acid_mutations", help="JSON of amino acid mutations per node and sequences for the root node")
    parser.add_argument("reconstructed_sequences", help="JSON with reconstructed sequences")

    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load nucleotide mutations.
    with open(args.nucleotide_mutations, "r") as fh:
        nuc_mutations = json.load(fh)

    # Load amino acid mutations.
    with open(args.amino_acid_mutations, "r") as fh:
        aa_mutations = json.load(fh)

    # Reconstruct sequences from mutations.
    reconstructed_sequences = {
        "annotations": aa_mutations["annotations"],
        "nodes": reconstruct_sequences_from_mutations(tree, nuc_mutations["nodes"], aa_mutations["nodes"])
    }

    # Export the reconstructed sequences to JSON.
    with open(args.reconstructed_sequences, "w") as oh:
        json.dump(reconstructed_sequences, oh, indent=2)
