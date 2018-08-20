"""
Plot the frequency trajectories of clades including:

  - observed frequencies over all time periods
  - censored frequencies at a specific timepoint
  - predicted frequencies from a specific timepoint and delta t

Inputs to these plots include:

  - a tree to estimate frequencies from
  - a set of fitness model parameters for a given set of predictors
"""
import argparse
import json
import os
import sys

# augur imports.
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "augur"))
from base.fitness_model import FitnessModel
from base.frequencies import KdeFrequencies

from .utils import load_tree_from_json_filename, load_frequencies_from_json_filename


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice JSON tree")
    parser.add_argument("frequencies", help="JSON with frequencies estimated from the given tree and used to estimate the given parameters")
    parser.add_argument("parameters", help="tab-delimited file of model parameters produced by the forecasting builds (e.g., model_parameters/.../0.tab)")

    args = parser.parse_args()

    # Load tree.
    tree = load_tree_from_json_filename(args.tree)

    # Load frequencies.
    frequencies = load_frequencies_from_json_filename(args.frequencies)

    # Populate the fitness model with the given tree, frequencies, and parameters.
    model = FitnessModel(
        tree,
        frequencies,
        time_interval,
        args.predictors,
        censor_frequencies=not args.no_censoring,
        epitope_masks_fname="%s/builds/flu/metadata/ha_masks.tsv" % code_directory,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        min_freq=0.1,
        predictor_kwargs=predictor_kwargs,
        sigma=args.sigma
    )
