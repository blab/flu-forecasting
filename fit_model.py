import argparse
import datetime
import json
import numpy as np
import os

# Add augur source to the Python path.
import sys
root_directory = os.path.dirname(os.path.abspath(__file__))
augur_directory = os.path.join(root_directory, "dist", "augur")
sys.path.insert(0, augur_directory)

# Load augur modules.
from base.fitness_model import fitness_model as FitnessModel, make_pivots
from base.frequencies import KdeFrequencies
from base.io_util import json_to_tree, json_to_clade_frequencies, reconstruct_sequences_from_tree_and_root
from base.process import process
from base.titer_model import TiterCollection


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
    parser.add_argument("--no-censoring", action="store_true", help="Disable censoring of future data during frequency estimation")
    parser.add_argument("--end-date", type=float, help="Maximum date to use data from when fitting the model")
    parser.add_argument("--step-size", type=float, default=0.5, help="Step size in years between timepoints the model fits to")
    parser.add_argument("--data-frame", help="optional name of a file to save the resulting model's data frame to")

    args = parser.parse_args()
    predictor_kwargs = {}

    # Load HA tree.
    ha_tree = load_tree_from_json_filename(args.ha_tree)

    # Load HA metadata.
    with open(args.ha_metadata, "r") as fh:
        ha_metadata = json.load(fh)

    # Get HA genes in order.
    ha_ordered_genes = get_ordered_genes_from_metadata(ha_metadata)

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

    frequencies = KdeFrequencies.from_json(json_frequencies)

    # Setup predictor arguments.
    predictor_kwargs = {
        "tau": 0.5,
        "time_window": 0.1
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
        censor_frequencies=not args.no_censoring,
        epitope_masks_fname="%s/builds/flu/metadata/ha_masks.tsv" % augur_directory,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        min_freq=0.1,
        predictor_kwargs=predictor_kwargs,
        end_date=args.end_date,
        step_size=args.step_size
    )
    model.predict()
    model.validate_prediction()

    # Save resulting model.
    model.to_json(args.model)

    # Save model's input data frame to a file, if a filename is given.
    if args.data_frame:
        df = model.to_data_frame()
        df.to_csv(args.data_frame, sep="\t", header=True, index=False)
