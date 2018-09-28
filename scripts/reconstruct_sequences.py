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


def reconstruct_sequences_from_mutations(tree, root_sequences):
    """Returns a tree for which each node is annotated with that node's
    corresponding nucleotide and amino acid sequences as reconstructed from the
    given root node's sequences and a tree with nucleotide and amino acid
    mutations annotated per node.
    """
    annotated_tree = copy.copy(tree)

    # Annotate root sequences.
    annotated_tree.root.sequences = root_sequences

    # Reconstruct sequences for all other nodes in the tree.
    for node in annotated_tree.find_clades():
        for child in node.clades:
            child_sequences = node.sequences.copy()

            # Merge mutations into a single data structure that can be iterated over once.
            mutation_sets = child.aa_muts.copy()
            mutation_sets["nuc"] = child.muts

            # Reconstruct amino acid sequences.
            for gene, mutations in mutation_sets.items():
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
            child.sequences = child_sequences

    return annotated_tree


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="tree JSON for auspice")
    parser.add_argument("sequences", help="root node's sequences JSON for auspice")
    parser.add_argument("annotated_tree", help="annotate tree JSON with reconstructed sequences")

    args = parser.parse_args()

    # Load tree.
    with open(args.tree, "r") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)

    # Load root sequences.
    with open(args.sequences, "r") as fh:
        root_sequences = json.load(fh)

    # Reconstruct sequences from mutations.
    annotated_tree = reconstruct_sequences_from_mutations(tree, root_sequences)

    # Export the annotated tree to JSON.
    extra_attr = [variable for variable in vars(annotated_tree).keys()
                  if variable not in ["clades", "up"]]
    annotated_json_tree = tree_to_json(annotated_tree, extra_attr=extra_attr)

    with open(args.annotated_tree, "w") as oh:
        json.dump(annotated_json_tree, oh, indent=2)
