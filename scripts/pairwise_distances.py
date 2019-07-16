"""Calculate the distance between sequences between seasons.
"""
import argparse

from augur.distance import read_distance_map, get_distance_between_nodes
from augur.frequency_estimators import TreeKdeFrequencies
from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, read_node_data, write_json

import Bio.Phylo
from collections import defaultdict
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--frequencies", help="frequencies JSON", required=True)
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", type=str, help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--attribute-name", nargs="+", help="name to store distances associated with the given distance map; multiple attribute names are linked to corresponding positional comparison method and distance map arguments", required=True)
    parser.add_argument("--map", nargs="+", help="JSON providing the distance map between sites and, optionally, sequences present at those sites; the distance map JSON minimally requires a 'default' field defining a default numeric distance and a 'map' field defining a dictionary of genes and one-based coordinates", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date", required=True)
    parser.add_argument("--years-back-to-compare", type=int, help="number of years prior to the current season to search for samples to calculate pairwise comparisons with", required=True)
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)

    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies_json = json.load(fh)

    frequencies = TreeKdeFrequencies.from_json(frequencies_json)
    pivots = frequencies.pivots

    # Identify pivots that belong within our search window for past samples.
    past_pivot_indices = (pivots < pivots[-1]) & (pivots >= pivots[-1] - args.years_back_to_compare)

    # Load sequences.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Index sequences by node name and gene.
    sequences_by_node_and_gene = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Load date annotations and annotate tree with them.
    date_annotations = read_node_data(args.date_annotations)
    for node in tree.find_clades():
        node.attr = date_annotations["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Identify samples to compare including those in the current timepoint
    # (pivot) and those in previous timepoints.
    current_samples = []
    past_samples = []
    date_by_sample = {}
    for tip in tree.find_clades(terminal=True):
        # Samples with nonzero frequencies in the last timepoint are current
        # samples. Those with one or more nonzero frequencies in the search
        # window of the past timepoints are past samples.
        if frequencies.frequencies[tip.name][-1] > 0:
            current_samples.append(tip.name)
        elif (frequencies.frequencies[tip.name][past_pivot_indices] > 0).sum() > 0:
            past_samples.append(tip.name)

        date_by_sample[tip.name] = tip.attr["numdate"]

    print("Expecting %i comparisons" % (len(current_samples) * len(past_samples) * len(args.attribute_name)))

    distances_by_node = {}
    distance_map_names = []
    comparisons = 0
    for attribute, distance_map_file in zip(args.attribute_name, args.map):
        # Load the given distance map.
        distance_map = read_distance_map(distance_map_file)
        distance_map_names.append(distance_map.get("name", distance_map_file))

        for current_sample in current_samples:
            if not current_sample in distances_by_node:
                distances_by_node[current_sample] = {}

            if not attribute in distances_by_node[current_sample]:
                distances_by_node[current_sample][attribute] = {}

            for past_sample in past_samples:
                # The past is in the past.
                comparisons += 1
                if date_by_sample[past_sample] < date_by_sample[current_sample]:
                    distances_by_node[current_sample][attribute][past_sample] = get_distance_between_nodes(
                        sequences_by_node_and_gene[past_sample],
                        sequences_by_node_and_gene[current_sample],
                        distance_map
                    )

    print("Calculated %i comparisons" % comparisons)
    # Prepare params for export.
    params = {
        "attribute": args.attribute_name,
        "map_name": distance_map_names,
        "years_back_to_compare": args.years_back_to_compare
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": distances_by_node}, args.output)
