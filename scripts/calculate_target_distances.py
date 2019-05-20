"""
Calculate pairwise distances between samples at adjacent timepoints (t and t - delta_t).
"""
import argparse
from augur.distance import read_distance_map, get_distance_between_nodes
from augur.frequency_estimators import TreeKdeFrequencies
from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, read_node_data, write_json
import Bio.Phylo
from collections import defaultdict
import json
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert distances JSON to a data frame",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--frequencies", help="frequencies JSON for all timepoints up to the current timepoint", required=True)
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", type=str, help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--attribute-name", help="name to store target distances", required=True)
    parser.add_argument("--map", help="JSON providing the distance map between sites and, optionally, sequences present at those sites; the distance map JSON minimally requires a 'default' field defining a default numeric distance and a 'map' field defining a dictionary of genes and one-based coordinates", required=True)
    parser.add_argument("--delta-pivots", type=int, default=2, help="number of frequency pivots to look back in time for pairwise distance calculations")
    parser.add_argument("--min-frequency", type=float, default=1e-5, help="minimum frequency of a sample to consider it for pairwise comparisons")
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)
    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    frequencies = TreeKdeFrequencies.from_json(frequencies_json)

    # Load sequences.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Index sequences by node name and gene.
    sequences_by_node_and_gene = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Load the given distance map.
    distance_map = read_distance_map(args.map)

    # Find all samples at non-zero frequencies in the current timepoint and the previous timepoint of interest.
    current_samples = []
    past_samples = []
    for tip in tree.find_clades(terminal=True):
        if frequencies.frequencies[tip.name][-1] > args.min_frequency:
            current_samples.append(tip.name)

        if frequencies.frequencies[tip.name][-(args.delta_pivots + 1)] > args.min_frequency:
            past_samples.append(tip.name)

    print(
        "Found %i current samples and %i past samples for %i total comparisons" % (len(current_samples), len(past_samples), (len(current_samples) * len(past_samples)))
    )

    # Calculate pairwise distances between all current nodes and also current
    # and past nodes.
    count = 0
    distances_by_node = {}
    for sample in current_samples:
        distances_by_node[sample] = {
            args.attribute_name: {}
        }

        for past_sample in past_samples + current_samples:
            distances_by_node[sample][args.attribute_name][past_sample] = get_distance_between_nodes(
                sequences_by_node_and_gene[past_sample],
                sequences_by_node_and_gene[sample],
                distance_map
            )
            count += 1

    # Write out the node annotations.
    print("Saving %i pairwise distances" % count)
    write_json({"nodes": distances_by_node}, args.output, indent=None)
