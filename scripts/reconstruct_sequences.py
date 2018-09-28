"""Calculate LBI for a given tree and one or more sets of parameters.
"""
import argparse
import copy
import json
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.io_util import json_to_tree, tree_to_json


def reconstruct_sequences_from_mutations(tree, nuc_mutations, aa_mutations):
    """Returns a dictionary of nucleotide and amino acid sequences per node as
    reconstructed from the given tree and mutations with the root node's
    sequences.

    Args:
        tree: a Bio.Phylo instance
        nuc_mutations: a dictionary of nucleotide sequences and mutations per node in the given tree
        aa_mutations: a dictionary of amino acid mutations per node with full length sequences annotated for the root node

    Returns:

        a dictionary of nucleotide and amino acid sequences per node in the
        given tree and the start/end coordinate annotations for the nucleotide
        and amino acid segments
    """
    # Annotate root sequences.
    sequences = {
        tree.root.name: {"nuc": nuc_mutations[tree.root.name]}
    }
    sequences[tree.root.name].update(aa_mutations[tree.root.name])

    # Reconstruct sequences for all other nodes in the tree.
    for node in tree.find_clades():
        for child in node.clades:
            # Copy the parent node's sequences as the default, assuming no
            # mutations have occurred.
            child_sequences = sequences[node.name].copy()

            # Annotate child's nucleotide sequences which already exist.
            child_sequences["nuc"] = nuc_mutations[child.name]

            # Reconstruct amino acid sequences.
            for gene, mutations in aa_mutations[child.name].items():
                if len(mutations) > 0:
                    # Convert sequence string to a list for in place manipulation.
                    gene_sequence = list(child_sequences[gene])

                    for mutation in mutations:
                        ancestral_aa = mutation[0]
                        derived_aa = mutation[-1]
                        position = int(mutation[1:-1])

                        assert gene_sequence[position - 1] == ancestral_aa
                        gene_sequence[position - 1] = derived_aa

                    # Convert list back to a string for the final child sequence.
                    child_sequences[gene] = "".join(gene_sequence)

                    assert child_sequences[gene] != node.sequences[gene]

            # Assign child sequences to child node.
            sequences[child.name] = child_sequences

    return sequences


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
        "nodes": reconstruct_sequences_from_mutations(tree, nuc_mutations, aa_mutations)
    }

    # Export the reconstructed sequences to JSON.
    with open(args.reconstructed_sequences, "w") as oh:
        json.dump(reconstructed_sequences, oh, indent=2)
