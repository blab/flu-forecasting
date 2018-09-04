"""Calculate the distance between amino acid sequences across entire genes or
at a predefined subset of sites.
"""
import argparse
import Bio
import Bio.Phylo
import json
import numpy as np
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.scores import mask_distance, read_masks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate amino acid sequence Hamming distance across entire genes or at a predefined subset of sites between root node of a tree and all other nodes in the tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("tree", help="Newick tree")
    parser.add_argument("translations", help="JSON of amino acid translations indexed by node name and gene name")
    parser.add_argument("output", help="JSON file with calculated distances stored by node name and attribute name")
    parser.add_argument("--attribute-names", nargs="+", help="names to store distances associated with the corresponding masks", required=True)
    parser.add_argument("--masks", help="tab-delimited mask definitions with mask name in first column and binary mask in second column")
    parser.add_argument("--mask-names", nargs="*", help="name of each mask to use from the given masks file. Distances are calculated for each given mask and stored in each of the corresponding attribute names. If no mask is provided, all sites will be used.")

    args = parser.parse_args()

    # Load tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Identify the name of the root node.
    root_node_name = tree.root.name

    # Load translation.
    with open(args.translations, "r") as json_fh:
        json_translations = json.load(json_fh)

    # Order gene annotations by their start position.
    # Annotations are dictionaries of start and end coordinates plus strand stored by gene name.
    annotations = sorted(
        json_translations["annotations"].items(),
        key=lambda item: item[1]["start"]
    )

    # Extract gene segments in order from annotations.
    # Only keep genes that have a corresponding translation in the root node.
    genes = [annotation[0]
             for annotation in annotations
             if annotation[0] in json_translations["nodes"][root_node_name]["translations"]]

    # Create a single amino acid sequence for the root node based on the order of the genes.
    # Convert the string into an array for mask comparisons.
    root_node_sequence = np.fromstring(
        "".join([
            json_translations["nodes"][root_node_name]["translations"][gene]
            for gene in genes
        ]),
        dtype="S1"
    )

    # Determine which sites to include in distance calculations.
    # If any masks are specified, calculate the distance for each mask and store it in the corresponding attribute name.
    if args.masks:
        # Load masks.
        masks = read_masks(args.masks)

        # Map masks to attribute names where distances will be stored.
        attributes_by_mask = dict(zip(args.mask_names, args.attribute_names))
    else:
        # If no masks are specified, calculate the distance using all sites.
        mask_name = "all_sites"
        masks = {mask_name: np.ones(len(root_node_sequence)).astype(bool)}

        # Use the requested attribute name to store distances for this mask.
        attributes_by_mask = {mask_name: args.attribute_names[0]}

    # Calculate Hamming distance between the root node's sequence and the single
    # amino acid sequence for each node in the tree.
    distances_by_node = {}
    for node, node_translations in json_translations["nodes"].items():
        node_sequence = np.fromstring(
            "".join([
                node_translations["translations"][gene]
                for gene in genes
            ]),
            dtype="S1"
        )

        # Calculate the distance for this node for each requested mask and store
        # the result in the corresponding attribute name.
        distances_by_node[node] = {}
        for mask_name, attribute in attributes_by_mask.items():
            distances_by_node[node][attribute] = mask_distance(
                root_node_sequence,
                node_sequence,
                masks[mask_name]
            )

    # Export distances to JSON.
    with open(args.output, "w") as oh:
        json.dump({"nodes": distances_by_node}, oh, indent=2)
