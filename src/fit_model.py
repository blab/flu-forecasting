import argparse
from augur.frequencies import TreeKdeFrequencies
from augur.titer_model import TiterCollection
from augur.utils import json_to_tree
from collections import OrderedDict
import datetime
import json
import numpy as np
import os
import sys

from forecast.fitness_model import fitness_model as FitnessModel, make_pivots, mean_absolute_error, sum_of_squared_errors


def reconstruct_sequences_from_tree_and_root(tree, root_sequences, ordered_genes):
    """Returns a tree for which each node is annotated with that node's
    corresponding nucleotide and amino acid sequences as reconstructed from the
    given root node's sequences and a tree with nucleotide and amino acid
    mutations annotated per node.

    The given sequence of gene names should be ordered by their genomic
    coordinates such that the annotated translations are stored in coordinate
    order.
    """
    # Annotate root translations using gene order information.
    tree.root.translations = OrderedDict()
    for gene in ordered_genes:
        tree.root.translations[gene] = root_sequences[gene]

    # Reconstruct sequences for all other nodes in the tree.
    for node in tree.find_clades():
        for child in node.clades:
            child_sequences = node.translations.copy()

            # Merge mutations into a single data structure that can be iterated over once.
            mutation_sets = child.aa_muts.copy()

            if "nuc" in ordered_genes:
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

    return tree


def load_tree_from_json_filename(filename):
    # Load JSON tree.
    with open(filename, "r") as json_fh:
         json_tree = json.load(json_fh)

    # Convert JSON tree layout to a Biopython Clade instance.
    tree = json_to_tree(json_tree)

    return tree


def pivot_to_date(pivot):
    """
    Convert the given pivot floating point date to its corresponding Python date instance.

    >>> pivot_to_date(2008.0)
    date(2008, 1, 1)
    >>> pivot_to_date(2008.75)
    date(2008, 10, 1)
    """
    year = int(pivot)
    month = int(round((pivot - year) * 12, 0)) + 1
    day = 1

    return datetime.date(year, month, day)


def get_ordered_genes_from_metadata(metadata_json):
    """Get gene names in order of their genomic coordinates. Annotations are stored
    in the metadata as a dictionary of start/end coordinates and strand indexed
    by gene name.
    """
    sorted_gene_annotations = sorted(
        metadata_json["annotations"].items(),
        key=lambda item: item[1]["start"]
    )
    ordered_genes = [gene for gene, annotations in sorted_gene_annotations]

    return ordered_genes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ha_tree", help="auspice tree JSON for HA")
    parser.add_argument("ha_metadata", help="auspice metadata JSON for HA")
    parser.add_argument("ha_sequences", help="auspice sequence JSON for HA")
    parser.add_argument("frequencies", help="JSON containing frequencies estimated from the given tree")
    parser.add_argument("model", help="output model JSON")
    parser.add_argument("predictors", nargs="+", help="one or more predictors to build model for")
    parser.add_argument("--na-tree", help="auspice tree JSON for NA")
    parser.add_argument("--na-sequences", help="auspice sequence JSON for NA")
    parser.add_argument("--titers", help="tab-delimited file of titer measurements")
    parser.add_argument("--dms", help="tab-delimited file of DMS preferences")
    parser.add_argument("--masks", help="tab-delimited file of mutational masks to use")
    parser.add_argument("--no-censoring", action="store_true", help="Disable censoring of future data during frequency estimation")
    parser.add_argument("--end-date", type=float, help="Maximum date to use data from when fitting the model")
    parser.add_argument("--step-size", type=float, default=0.5, help="Step size in years between timepoints the model fits to")
    parser.add_argument("--delta-time", type=float, default=1.0, help="Number of years to project forward from each timepoint")
    parser.add_argument("--tip-data-frame", help="optional name of a file to save the resulting model's tip data frame to")
    parser.add_argument("--clade-data-frame", help="optional name of a file to save the resulting model's clade data frame to")
    parser.add_argument("--validation-data-frame", help="optional name of a file to save the resulting model's cross-validation results to")
    parser.add_argument("--prepare-only", action="store_true", help="prepare model inputs without fitting model parameters")
    parser.add_argument("--min-freq", type=float, default=0.1, help="minimum frequency for clades to be used in model fitting")
    parser.add_argument("--max-freq", type=float, default=0.99, help="maximum frequency for clades to be used in model fitting")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()
    predictor_kwargs = {}

    # Load HA tree.
    ha_tree = load_tree_from_json_filename(args.ha_tree)

    # Load HA metadata.
    with open(args.ha_metadata, "r") as fh:
        ha_metadata = json.load(fh)

    # Get HA genes in order excluding nucleotide sequences.
    ha_ordered_genes = [
        gene
        for gene in get_ordered_genes_from_metadata(ha_metadata)
        if gene != "nuc"
    ]

    # Load HA root sequences.
    with open(args.ha_sequences, "r") as fh:
        ha_root_sequences = json.load(fh)

    # Reconstruct ordered translations from the tree, root sequences, and gene order.
    ha_tree = reconstruct_sequences_from_tree_and_root(ha_tree, ha_root_sequences, ha_ordered_genes)

    # Load NA tree if it has been provided.
    if args.na_tree:
        na_tree = load_tree_from_json_filename(args.na_tree)

        # Annotate epitope mutations from NA tree onto HA tree.
        na_mutations_by_strain = {
            node.name: node.attr["ep"]
            for node in na_tree.find_clades()
            if node.is_terminal()
        }

        for node in ha_tree.find_clades():
            if node.is_terminal():
                node.attr["na_ep"] = na_mutations_by_strain.get(node.name, 0)

    # Load frequencies.
    with open(args.frequencies, "r") as json_fh:
        json_frequencies = json.load(json_fh)

    frequencies = TreeKdeFrequencies.from_json(json_frequencies)

    # Setup predictor arguments.
    predictor_kwargs = {
        "tau": 0.5,
        "time_window": 0.1,
        "step_size": args.step_size
    }

    # If DMS preferences were provided, link to them.
    if args.dms:
        predictor_kwargs["preferences_file"] = args.dms

    # If titers were provided, load them for the model to use.
    if args.titers:
        predictor_kwargs.update({
            "lam_avi": 2.0,
            "lam_pot": 0.3,
            "lam_drop": 2.0
        })
        titers, strains, sources = TiterCollection.load_from_file(args.titers)
        predictor_kwargs["titers"] = titers

    # Run the fitness model.
    model = FitnessModel(
        ha_tree,
        frequencies,
        args.predictors,
        cross_validate=True,
        censor_frequencies=not args.no_censoring,
        epitope_masks_fname=args.masks,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        min_freq=args.min_freq,
        max_freq=args.max_freq,
        predictor_kwargs=predictor_kwargs,
        end_date=args.end_date,
        step_size=args.step_size,
        delta_time=args.delta_time,
        verbose=int(args.verbose),
        enforce_positive_predictors=False,
        cost_function=sum_of_squared_errors
    )

    if args.prepare_only:
        model.prep_nodes()
        model.calc_node_frequencies()
        model.calc_all_predictors()
        model.standardize_predictors()
        model.select_nonoverlapping_clades_for_fitting()

        with open(args.model, "w") as fh:
            fh.write("{}\n")
    else:
        validation_df = model.predict()
        model.validate_prediction()
        model.validate_prediction(test=True)

        # Save resulting model.
        model.to_json(args.model)

        # Save the model's cross-validation results.
        validation_df.to_csv(args.validation_data_frame, sep="\t", index=False)

    # Save model's input data frame to a file, if a filename is given.
    if args.tip_data_frame:
        tip_df, clade_df = model.to_data_frames()
        tip_df.to_csv(args.tip_data_frame, sep="\t", header=True, index=False)
        clade_df.to_csv(args.clade_data_frame, sep="\t", header=True, index=False)
