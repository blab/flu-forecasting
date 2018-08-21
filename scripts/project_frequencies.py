"""
Project the frequency trajectories of clades including:

  - observed frequencies over all time periods
  - censored frequencies at a specific timepoint
  - predicted frequencies from a specific timepoint and delta t

Inputs to these plots include:

  - a tree to estimate frequencies from
  - a set of fitness model parameters for a given set of predictors
"""
import argparse
from collections import defaultdict
import json
import numpy as np
import os
import pandas as pd
import sys

import matplotlib.pyplot as plt

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.fitness_model import fitness_model as FitnessModel
from base.frequencies import KdeFrequencies

from utils import load_tree_from_json_filename, load_frequencies_from_json_filename


def project_clade_frequencies_by_delta_from_time(tree, model, time, delta):
    """
    Project clade frequencies from a given time to the future by a given delta.
    """
    # Calculate the steps between the projection date and delta time into the
    # future. First, find the frequency pivot that is closest to the requested
    # projection date.
    max_date = model.timepoints[np.searchsorted(model.timepoints, time)]
    future_date = max_date + delta

    # Then, calculate a fixed number of steps between that pivot and delta time
    # into the future.
    projected_pivots = np.linspace(max_date, future_date, int(delta_steps_per_year * delta))
    deltas = projected_pivots - max_date

    # Identify tip predictors and frequencies at the current time point.
    all_pred = model.predictor_arrays[max_date]
    all_freqs = model.freq_arrays[max_date]

    # For each requested delta, project current tip frequencies using the model
    # and calculate the corresponding projected clade frequencies.
    projected_clade_frequencies = defaultdict(list)

    for delta in deltas:
        # Project all tip frequencies.
        pred_freq = model.projection(model.model_params, all_pred, all_freqs, delta)

        # Normalize projected frequencies.
        pred_freq = pred_freq / pred_freq.sum()

        # Store projected frequencies by clade id.
        for i, tip in enumerate(model.tips):
            projected_clade_frequencies[tip.clade].append(pred_freq[i])

        # Calculate projected frequencies for internal nodes and store by clade it.
        for node in tree.find_clades(order="postorder"):
            if not node.is_terminal():
                projected_clade_frequencies[node.clade].append(pred_freq[node.tips].sum())

    projected_frequencies = {
        "params": {
            "max_date": max_date
        },
        "data": {
            "pivots": projected_pivots.tolist(),
            "frequencies": projected_clade_frequencies
        }
    }

    return projected_frequencies


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("tree", help="auspice JSON tree")
    parser.add_argument("frequencies", help="JSON with frequencies estimated from the given tree and used to estimate the given parameters")
    parser.add_argument("parameters", help="tab-delimited file of model parameters produced by the forecasting builds (e.g., model_parameters/.../0.tab)")
    parser.add_argument("projected_frequencies", help="JSON with frequencies projected into the future")
    parser.add_argument("projection_date", type=float, help="date from which frequencies should be projected")
    parser.add_argument("delta", type=float, default=1.0, help="amount of time in years to project frequencies into the future")
    delta_steps_per_year = 12

    args = parser.parse_args()

    # Load tree.
    tree = load_tree_from_json_filename(args.tree)

    # Load frequencies.
    frequencies = load_frequencies_from_json_filename(args.frequencies)

    # Load model parameters.
    parameters_df = pd.read_table(args.parameters)
    predictors = {record["predictor"]: (record["param"], record["global_sd"])
                  for record in parameters_df.to_dict("records")}

    # Setup predictor arguments.
    predictor_kwargs = {}

    # Populate the fitness model with the given tree, frequencies, and parameters.
    model = FitnessModel(
        tree,
        frequencies,
        predictors,
        epitope_masks_fname="%s/builds/flu/metadata/ha_masks.tsv" % augur_path,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        min_freq=0.1,
        predictor_kwargs=predictor_kwargs
    )
    model.predict()

    # Calculate projected frequencies.
    projected_frequencies = project_clade_frequencies_by_delta_from_time(
        tree,
        model,
        args.projection_date,
        args.delta
    )

    # Export projected frequencies in the context of the input frequencies.
    with open(args.projected_frequencies, "w") as oh:
        json_frequencies = json.dump(projected_frequencies, oh, indent=2)

    # Calculate censored frequencies.
    frequency_parameters = frequencies.get_params()
    frequency_parameters["max_date"] = projected_frequencies["params"]["max_date"]
    censored_frequencies = KdeFrequencies(**frequency_parameters)
    censored_frequencies.estimate(tree)

    #
    # Plot frequencies
    #

    # Pick a clade to plot.
    clade = [node for node in tree.get_nonterminals()
             if len(node.get_terminals()) > 20 and len(node.get_terminals()) < 100][2]
    clade_id = clade.clade

    # Plot observed frequencies.
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(frequencies.pivots, frequencies.frequencies[clade_id], "o-", label="Observed")

    # Plot censored frequencies.
    ax.plot(censored_frequencies.pivots, censored_frequencies.frequencies[clade_id], "o-", label="Censored")

    # Plot projected frequencies.
    ax.plot(projected_frequencies["data"]["pivots"], projected_frequencies["data"]["frequencies"][clade_id], "o--", label="Projected")

    # Annotate projection date.
    ax.axvline(x=projected_frequencies["params"]["max_date"], color="#000000", alpha=0.5, zorder=-1)

    # Label axes.
    ax.set_xlabel("Date")
    ax.set_ylabel("Frequency")

    # Annotate clade membership and size.
    ax.set_title("Clade from %s with %i tips" % (clade.attr["clade_membership"], len(clade.get_terminals())))

    ax.legend()

    plt.savefig("analyses/projected-frequencies/example-clade.png")
    plt.close()
