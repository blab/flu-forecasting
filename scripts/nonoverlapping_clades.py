"""Find clades in a given tree by distinct haplotypes in the given amino acid sequences corresponding to internal nodes in the tree.
"""
import argparse
from augur.frequency_estimators import TreeKdeFrequencies
from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, write_json
import Bio.Phylo
import Bio.SeqIO
import hashlib
import json
import pandas as pd

# Magic number of maximum length of SHA hash to keep for each clade.
MAX_HASH_LENGTH = 7


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find clades in a tree by distinct amino acid haplotypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick tree to identify clades in")
    parser.add_argument("--translations", required=True, nargs="+", help="FASTA file(s) of amino acid sequences per node")
    parser.add_argument("--gene-names", required=True, nargs="+", help="gene names corresponding to translations provided")
    parser.add_argument("--output", required=True, help="JSON of clade annotations for nodes in the given tree")
    parser.add_argument("--output-tip-clade-table", help="optional table of all clades per tip in the tree")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the tip clade output table in the format of 'key=value' pairs")

    args = parser.parse_args()

    # Load the tree.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load translations for nodes in the given tree and index them by gene name and node name.
    translations = load_alignments(args.translations, args.gene_names)
    translations_by_gene_name = {}
    for gene in translations:
        translations_by_gene_name[gene] = {}
        for seq in translations[gene]:
            translations_by_gene_name[gene][seq.name] = str(seq.seq)

    clades = {}
    for node in tree.find_clades(order="preorder", terminal=False):
        # Assign the current node a clade id based on the hash of its
        # full-length amino acid sequence.
        node_sequence = "".join([translations_by_gene_name[gene][node.name] for gene in args.gene_names])
        clades[node.name] = {"clade_membership": hashlib.sha256(node_sequence.encode()).hexdigest()[:MAX_HASH_LENGTH]}

        # Assign the current node's clade id to all of its terminal children.
        for child in node.clades:
            if child.is_terminal():
                clades[child.name] = clades[node.name]

    # Count unique clade groups.
    distinct_clades = {clade["clade_membership"] for clade in clades.values()}
    print("Found %i distinct clades" % len(distinct_clades))

    # Write out the node annotations.
    write_json({"nodes": clades}, args.output)

    # Output the optional tip-to-clade table, if requested.
    if args.output_tip_clade_table:
        records = []
        for tip in tree.find_clades(terminal=True):
            # Note the tip's own clade assignment which may be distinct from its
            # parent's.
            depth = 0
            records.append([tip.name, clades[tip.name]["clade_membership"], depth])

            parent = tip.parent
            depth += 1
            while True:
                records.append([tip.name, clades[parent.name]["clade_membership"], depth])

                if parent == tree.root:
                    break

                parent = parent.parent
                depth += 1

        df = pd.DataFrame(records, columns=["tip", "clade_membership", "depth"])
        df = df.drop_duplicates(subset=["tip", "clade_membership"])

        # Add any additional annotations requested by the user in the format of
        # "key=value" pairs where each key becomes a new column with the given
        # value.
        if args.annotations:
            for annotation in args.annotations:
                key, value = annotation.split("=")
                df[key] = value

        df.to_csv(args.output_tip_clade_table, sep="\t", index=False)
