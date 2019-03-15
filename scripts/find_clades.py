"""Find clades in a given tree by distinct haplotypes in the given amino acid sequences corresponding to internal nodes in the tree.
"""
import argparse
from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, write_json
import Bio.Phylo
import hashlib
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
    parser.add_argument("--use-incremental-ids", action="store_true", help="report an incremental integer id for each clade instead of concatenated mutations")
    parser.add_argument("--use-hash-ids", action="store_true", help="report an abbreviated SHA1 hash for each clade instead of concatenated mutations")
    parser.add_argument("--minimum-tips", type=int, default=1, help="minimum number of tips required for a clade to be assigned its own annotation")
    parser.add_argument("--output", required=True, help="JSON of clade annotations for nodes in the given tree")
    parser.add_argument("--output-tip-clade-table", help="optional table of all clades per tip in the tree")

    args = parser.parse_args()

    # Load the tree.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load translations and index them by gene name and node name.
    translations = load_alignments(args.translations, args.gene_names)
    translations_by_gene_name = {}
    for gene in translations:
        translations_by_gene_name[gene] = {}
        for seq in translations[gene]:
            translations_by_gene_name[gene][seq.name] = str(seq.seq)

    # Annotate mutations between each node and the MRCA of all nodes.
    clades = {}
    for node in tree.find_clades():
        if node == tree.root:
            clades[node.name] = {"clade_membership": "root"}
        elif node.count_terminals() < args.minimum_tips:
            # Assign tips and small clades to the same clade annotation as their
            # immediate parent.
            clades[node.name] = clades[node.parent.name].copy()
        else:
            # Collect all mutations between the current node and the MRCA.
            mutations = []
            for gene in args.gene_names:
                for i in range(len(translations_by_gene_name[gene][tree.root.name])):
                    if translations_by_gene_name[gene][tree.root.name][i] != translations_by_gene_name[gene][node.name][i]:
                        # Store mutations with gene, position, and allele like "HA1:131K".
                        mutations.append(f"{gene}:{i + 1}{translations_by_gene_name[gene][node.name][i]}")

            if len(mutations) > 0:
                # If this node has any mutations, concatenate them into a
                # comma-delimited string.
                clades[node.name] = {"clade_membership": ",".join(mutations)}
            else:
                # Otherwise, use the clade annotation of this node's parent.
                clades[node.name] = clades[node.parent.name].copy()

    # TODO: Remove clade annotations for internal nodes that have no tips
    # sharing the same annotation?

    # Count unique clade groups.
    distinct_clades = {clade["clade_membership"] for clade in clades.values()}
    print("Found %i distinct clades" % len(distinct_clades))

    # Replace concatenated mutations with incremental ids, if requested.
    if args.use_incremental_ids:
        # Assign clade numbers from root to tip.
        clade_number = 0
        mutations_to_number = {}
        for node in tree.find_clades():
            if clades[node.name]["clade_membership"] not in mutations_to_number:
                mutations_to_number[clades[node.name]["clade_membership"]] = clade_number
                clade_number += 1

        for node in tree.find_clades():
            clades[node.name]["clade_membership"] = "Clade %i" % mutations_to_number[clades[node.name]["clade_membership"]]
    elif args.use_hash_ids:
        # Assign abbreviated SHA hashes based on concatenated mutations.
        for node_name in clades.keys():
            if clades[node_name]["clade_membership"] != "root":
                clades[node_name]["clade_membership"] = hashlib.sha256(clades[node_name]["clade_membership"].encode()).hexdigest()[:MAX_HASH_LENGTH]

    # Write out the node annotations.
    write_json({"nodes": clades}, args.output)

    # Output the optional tip-to-clade table, if requested.
    if args.output_tip_clade_table:
        records = []
        for tip in tree.find_clades(terminal=True):
            parent = tip.parent
            depth = 1
            while True:
                records.append([tip.name, clades[parent.name]["clade_membership"], depth])

                if parent == tree.root:
                    break

                parent = parent.parent
                depth += 1

        df = pd.DataFrame(records, columns=["tip", "clade_membership", "depth"])
        df = df.drop_duplicates(subset=["tip", "clade_membership"])
        df.to_csv(args.output_tip_clade_table, sep="\t", index=False)
