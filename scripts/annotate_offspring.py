"""Annotate the standardized observed number of offspring per node to tips as a representation of observed fitness.
"""
import argparse
from augur.utils import annotate_parents_for_tree
import Bio.Phylo
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate the standardized observed number of offspring per node to tips as a representation of observed fitness",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", required=True, help="Newick file for the tree used to construct the given node data JSONs")
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints")
    parser.add_argument("--output", required=True, help="table of standardized offspring per sample")

    args = parser.parse_args()

    # Load the final tree.
    tree = Bio.Phylo.read(args.tree, "newick")

    # Annotate offspring per internal node.
    tree = annotate_parents_for_tree(tree)

    # Count offspring per node.
    offspring = []
    for node in tree.find_clades():
        if node.parent:
            parent_count = node.parent.offspring
        else:
            parent_count = 0

        node.offspring = parent_count + node.count_terminals()

        if node.is_terminal():
            offspring.append({
                "strain": node.name,
                "offspring": node.offspring
            })

    offspring_df = pd.DataFrame(offspring)

    # Load tip attributes per timepoint.
    tips = pd.read_csv(args.tip_attributes, sep="\t")

    # Annotate tips with offspring.
    tips = tips.merge(offspring_df, on="strain")

    # Calculate statistics for offspring per timepoint.
    mean_and_std_offspring = tips.groupby("timepoint").aggregate({
        "offspring": ["mean", "std"]
    }).reset_index().rename(columns={
        "mean": "mean_offspring",
        "std": "std_offspring"
    })
    mean_and_std_offspring.columns = ["timepoint", "mean_offspring", "std_offspring"]

    # Annotate tips with mean and std offspring.
    tips = tips.merge(mean_and_std_offspring, on="timepoint")

    # Standardize offspring per timepoint.
    tips["offspring"] = (tips["offspring"] - tips["mean_offspring"]) / tips["std_offspring"]

    # Drop intermediate columns.
    tips = tips.drop(columns=["mean_offspring", "std_offspring"]).copy()

    # Save tips with observed offspring annotations.
    tips.to_csv(args.output, sep="\t", header=True, index=False)
