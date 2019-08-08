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
    parser.add_argument("--frequencies", required=True, help="frequencies JSON for the given tree")
    parser.add_argument("--reference", required=True, help="GenBank file of reference sample to identify mutations against; should contain gene annotations to enable translation of nucleotide sequences to amino acids")
    parser.add_argument("--translations", required=True, nargs="+", help="FASTA file(s) of amino acid sequences per node")
    parser.add_argument("--gene-names", required=True, nargs="+", help="gene names corresponding to translations provided")
    parser.add_argument("--use-incremental-ids", action="store_true", help="report an incremental integer id for each clade instead of concatenated mutations")
    parser.add_argument("--use-hash-ids", action="store_true", help="report an abbreviated SHA1 hash for each clade instead of concatenated mutations")
    parser.add_argument("--minimum-tips", type=int, default=1, help="minimum number of tips required for a clade to be assigned its own annotation")
    parser.add_argument("--minimum-frequency", type=float, default=0.01, help="minimum frequency for a clade to be considered")
    parser.add_argument("--output", required=True, help="JSON of clade annotations for nodes in the given tree")
    parser.add_argument("--output-tip-clade-table", help="optional table of all clades per tip in the tree")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the tip clade output table in the format of 'key=value' pairs")

    args = parser.parse_args()

    # Load the tree.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    kde_frequencies = TreeKdeFrequencies.from_json(frequencies_json)

    # Load the reference and translate its sequences.
    reference = Bio.SeqIO.read(args.reference, "genbank")
    cds_features = [
        feature
        for feature in reference.features
        if feature.type == "CDS" and "gene" in feature.qualifiers
    ]

    reference_sequences_by_gene = {}
    for feature in cds_features:
        reference_sequences_by_gene[feature.qualifiers["gene"][0]] = str(feature.translate(reference).seq)

    # Load translations for nodes in the given tree and index them by gene name and node name.
    translations = load_alignments(args.translations, args.gene_names)
    translations_by_gene_name = {}
    for gene in translations:
        translations_by_gene_name[gene] = {}
        for seq in translations[gene]:
            translations_by_gene_name[gene][seq.name] = str(seq.seq)

    # Annotate mutations between each node and the given reference sample.
    clades = {}
    for node in tree.find_clades():
        if node == tree.root:
            clades[node.name] = {"clade_membership": "root"}
        elif node.count_terminals() < args.minimum_tips or kde_frequencies.frequencies[node.name][-1] < args.minimum_frequency:
            # Assign tips and small clades to the same clade annotation as their
            # immediate parent.
            clades[node.name] = clades[node.parent.name].copy()
        else:
            # If using hash ids, assign clades based on hashes of the
            # full-length amino acid sequence for each node. Otherwise, use
            # clades based on mutations relative to the given reference.
            if args.use_hash_ids:
                node_sequence = "".join([translations_by_gene_name[gene][node.name] for gene in args.gene_names])
                clades[node.name] = {"clade_membership": hashlib.sha256(node_sequence.encode()).hexdigest()[:MAX_HASH_LENGTH]}
            else:
                # Collect all mutations between the current node and the MRCA.
                mutations = []
                for gene in args.gene_names:
                    for i in range(len(reference_sequences_by_gene[gene])):
                        if reference_sequences_by_gene[gene][i] != translations_by_gene_name[gene][node.name][i]:
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
    # elif args.use_hash_ids:
    #     # Assign abbreviated SHA hashes based on concatenated mutations.
    #     for node_name in clades.keys():
    #         if clades[node_name]["clade_membership"] != "root":
    #             clades[node_name]["clade_membership"] = hashlib.sha256(clades[node_name]["clade_membership"].encode()).hexdigest()[:MAX_HASH_LENGTH]

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
