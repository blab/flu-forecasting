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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.fitness_model import fitness_model as FitnessModel
from base.frequencies import KdeFrequencies

from utils import load_tree_from_json_filename, load_frequencies_from_json_filename


def project_clade_frequencies_by_delta_from_time(tree, model, time, delta, delta_steps_per_year=12):
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
    parser.add_argument("model", help="fitness model JSON with learned parameters, projected frequencies, and predictor configuration (e.g., models/.../0.json)")
    parser.add_argument("projected_frequencies", help="JSON with frequencies projected into the future")
    parser.add_argument("delta", type=float, default=1.0, help="amount of time in years to project frequencies into the future")
    parser.add_argument("projection_dates", type=float, nargs="+", help="dates from which frequencies should be projected")

    args = parser.parse_args()

    # Load tree.
    tree = load_tree_from_json_filename(args.tree)

    # Load frequencies.
    frequencies = load_frequencies_from_json_filename(args.frequencies)
    frequencies.include_internal_nodes = True
    frequencies.estimate(tree)

    # Load the model.
    with open(args.model, "r") as fh:
        json_model = json.load(fh)

    predictors = {record["predictor"]: [round(record["param"], 2), round(record["global_sd"], 2)]
                  for record in json_model["params"]}
    predictors_key = "-".join(sorted([record["predictor"] for record in json_model["params"]]))

    # Setup predictor arguments.
    predictor_kwargs = json_model["predictor_kwargs"]

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
    model.prep_nodes()
    model.delta_time = json_model["delta_time"]

    predictor_arrays = {}
    for key in json_model["predictor_arrays"]:
        predictor_arrays[float(key)] = np.array(json_model["predictor_arrays"][key])

    model.predictor_arrays = predictor_arrays

    freq_arrays = {}
    for key in json_model["freq_arrays"]:
        freq_arrays[float(key)] = np.array(json_model["freq_arrays"][key])

    model.freq_arrays = freq_arrays

    model.select_clades_for_fitting()

    # Calculate projected frequencies for each requested date.
    records = []
    for projection_date in args.projection_dates:
        print("Project from %s" % projection_date)
        projected_frequencies = project_clade_frequencies_by_delta_from_time(
            tree,
            model,
            projection_date,
            args.delta
        )
        for node in model.tree.find_clades():
            for pivot, frequency in zip(projected_frequencies["data"]["pivots"],
                                        projected_frequencies["data"]["frequencies"][node.clade]):
                records.append({
                    "predictors": predictors_key,
                    "clade": node.clade,
                    "clade_membership": node.attr["clade_membership"],
                    "projection_date": projected_frequencies["params"]["max_date"],
                    "pivot": pivot,
                    "frequency": frequency
                })

    projected_df = pd.DataFrame(records)
    projected_df.to_csv(args.projected_frequencies, sep="\t", header=True, index=False)

    # # Export projected frequencies in the context of the input frequencies.
    # with open(args.projected_frequencies, "w") as oh:
    #     json_frequencies = json.dump(projected_frequencies, oh, indent=2)

    #
    # Plot frequencies
    #

    projection_dates = projected_df["projection_date"].unique()
    fig, axes = plt.subplots(
        len(projection_dates),
        1,
        figsize=(8, 2.5 * len(projection_dates)),
        squeeze=False,
        sharex=True,
        sharey=True,
        gridspec_kw={"hspace": 0.2}
    )
    print(projection_dates)

    for i, projection_date in enumerate(projection_dates):
        ax = axes.flatten()[i]

        # Pick a clade to plot.
        clade = model.fit_clades[projection_date][0]
        clade_id = clade.clade

        # Annotate clade membership and size.
        ax.set_title(
            "Clade from %s with %i tips (%s)" %
            (clade.attr["clade_membership"], len(clade.get_terminals()), predictors)
        )

        df = projected_df[(projected_df["projection_date"] == projection_date) &
                          (projected_df["clade"] == clade_id)]

        # Calculate mean +/- std of fitness for this clade vs. others.
        clade_fitness = model.predictor_arrays[projection_date][clade.tips]
        non_clade_tips = np.array([i for i in range(len(model.tips)) if i not in clade.tips])
        non_clade_fitness = model.predictor_arrays[projection_date][non_clade_tips]
        ax.text(0.05, 0.9, "Clade fitness: %.2f +/- %.2f" % (clade_fitness.mean(), clade_fitness.std()), transform=ax.transAxes)
        ax.text(0.05, 0.8, "Non-clade fitness: %.2f +/- %.2f" % (non_clade_fitness.mean(), non_clade_fitness.std()), transform=ax.transAxes)

        # Calculate censored frequencies.
        frequency_parameters = frequencies.get_params()
        frequency_parameters["max_date"] = projection_date
        censored_frequencies = KdeFrequencies(**frequency_parameters)
        censored_frequencies.estimate(tree)

        # Plot observed frequencies.
        ax.plot(frequencies.pivots, frequencies.frequencies[clade_id], "o-", label="Observed")

        # Plot censored frequencies.
        ax.plot(censored_frequencies.pivots, censored_frequencies.frequencies[clade_id], "o-", label="Censored")

        # Plot projected frequencies.
        ax.plot(df["pivot"], df["frequency"], "o--", label="Projected")

        # Annotate projection date.
        ax.axvline(x=projection_date, color="#000000", alpha=0.5, zorder=-1)

        # Label axes.
        ax.set_xlabel("Date")
        ax.set_ylabel("Frequency")

    ax.legend()

    plt.savefig("analyses/projected-frequencies/example-clade.png")
    plt.close()
