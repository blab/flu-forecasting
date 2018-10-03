"""Calculate LBI for a given tree and one or more sets of parameters.
"""
import argparse
import Bio
import Bio.Phylo
from collections import OrderedDict
import copy
import json
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.io_util import json_to_tree, tree_to_json


def reconstruct_sequences_from_tree_and_root(tree, root_sequences, ordered_genes):
    """Returns a tree for which each node is annotated with that node's
    corresponding nucleotide and amino acid sequences as reconstructed from the
    given root node's sequences and a tree with nucleotide and amino acid
    mutations annotated per node.

    The given sequence of gene names should be ordered by their genomic
    coordinates such that the annotated translations are stored in coordinate
    order.
    """
    annotated_tree = copy.deepcopy(tree)

    # Annotate root translations using gene order information.
    annotated_tree.root.translations = OrderedDict()
    for gene in ordered_genes:
        annotated_tree.root.translations[gene] = root_sequences[gene]

    # Reconstruct sequences for all other nodes in the tree.
    for node in annotated_tree.find_clades():
        for child in node.clades:
            child_sequences = node.translations.copy()

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

                    assert child_sequences[gene] != node.translations[gene]

            # Assign child sequences to child node.
            child.translations = child_sequences

    return annotated_tree


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
        tree.root.name: {"translations": {"nuc": nuc_mutations[tree.root.name]["sequence"]}}
    }
    sequences[tree.root.name]["translations"].update(aa_mutations[tree.root.name]["aa_sequences"])

    # Reconstruct sequences for all other nodes in the tree.
    for node in tree.find_clades():
        for child in node.clades:
            # Copy the parent node's sequences as the default, assuming no
            # mutations have occurred.
            child_sequences = copy.deepcopy(sequences[node.name])

            # Annotate child's nucleotide sequences which already exist.
            child_sequences["translations"]["nuc"] = nuc_mutations[child.name]["sequence"]

            # Reconstruct amino acid sequences.
            for gene, mutations in aa_mutations[child.name]["aa_muts"].items():
                if len(mutations) > 0:
                    # Convert sequence string to a list for in place manipulation.
                    gene_sequence = list(child_sequences["translations"][gene])

                    for mutation in mutations:
                        ancestral_aa = mutation[0]
                        derived_aa = mutation[-1]
                        position = int(mutation[1:-1])

                        assert gene_sequence[position - 1] == ancestral_aa
                        gene_sequence[position - 1] = derived_aa

                    # Convert list back to a string for the final child sequence.
                    child_sequences["translations"][gene] = "".join(gene_sequence)

                    assert child_sequences["translations"][gene] != sequences[node.name]["translations"][gene]

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
        "nodes": reconstruct_sequences_from_mutations(tree, nuc_mutations["nodes"], aa_mutations["nodes"])
    }

    # Export the reconstructed sequences to JSON.
    with open(args.reconstructed_sequences, "w") as oh:
        json.dump(reconstructed_sequences, oh, indent=2)
