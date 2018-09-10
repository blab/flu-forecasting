"""Calculate LBI for a given tree and one or more sets of parameters.
"""
import argparse
import Bio
import Bio.Phylo
from collections import defaultdict
import json
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.scores import select_nodes_in_season, calculate_LBI


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate LBI for a given tree and one or more sets of parameters.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("tree", help="Newick tree")
    parser.add_argument("branch_lengths", help="JSON with branch lengths and internal node dates estimated by TreeTime")
    parser.add_argument("output", help="JSON file with calculated distances stored by node name and attribute name")
    parser.add_argument("--attribute-names", nargs="+", help="names to store distances associated with the corresponding masks", required=True)
    parser.add_argument("--tau", nargs="+", type=float, help="tau value(s) defining the neighborhood of each clade", required=True)
    parser.add_argument("--window", nargs="+", type=float, help="time window(s) to calculate LBI across", required=True)

    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Load branch lengths.
    with open(args.branch_lengths, "r") as json_fh:
        branch_lengths = json.load(json_fh)

    # Annotate branch lengths and dates onto tree nodes.
    for node in tree.find_clades():
        node.attr = branch_lengths["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Find maximum time point in the given tree.
    timepoint = max(node.attr["num_date"] for node in tree.find_clades())

    # Calculate LBI for all requested sets of parameters and annotate LBI values
    # to the corresponding attribute names.
    lbi_by_node = defaultdict(dict)
    for i in range(len(args.tau)):
        tau = args.tau[i]
        window = args.window[i]
        attribute_name = args.attribute_names[i]

        # Select nodes that are alive in the given time window.
        select_nodes_in_season(tree, timepoint, window)

        # Calculate LBI.
        calculate_LBI(tree, attribute_name, tau)

        # Collect LBI values into a per-node JSON for export.
        for node in tree.find_clades():
            lbi_by_node[node.name][attribute_name] = node.attr[attribute_name]

    # Export LBI to JSON.
    with open(args.output, "w") as oh:
        json.dump({"nodes": lbi_by_node}, oh, indent=2)
