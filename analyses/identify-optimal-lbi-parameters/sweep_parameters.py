"""
Identify optimal LBI parameters

Perform a parameter sweep across a range of values for both the tau and time window parameters, fitting a fitness model for each parameter set, and identify the set that optimizes the fitness model's frequency correlation metric.
"""
import argparse

# Configure environment

import os
import sys

# Project path
project_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Path to top-level of the augur source code repository.
augur_path = os.path.join(project_path, "dist", "augur")
print(augur_path)
sys.path.append(augur_path)

import Bio
import Bio.Phylo
from collections import defaultdict
from datetime import datetime, timedelta
import itertools
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns

from base.io_util import json_to_tree
from base.frequencies import KdeFrequencies
from base.fitness_model import fitness_model as FitnessModel, get_matthews_correlation_coefficient_for_data_frame
from base.titer_model import TiterCollection

# Configure matplotlib theme.
fontsize = 14
matplotlib_params = {
    'axes.labelsize': fontsize,
    'font.size': fontsize,
    'legend.fontsize': 12,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'text.usetex': False,
    'figure.figsize': [8, 6],
    'savefig.dpi': 300,
    'figure.dpi': 300,
    'text.usetex': False
}

plt.rcParams.update(matplotlib_params)

# Turn off spines for all plots.
plt.rc("axes.spines", top=False, right=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice JSON tree")
    parser.add_argument("frequencies", help="JSON with frequencies estimated for the given tree")
    parser.add_argument("results", help="tab-delimited model results for all LBI parameters")

    args = parser.parse_args()

    # Load tree
    with open(args.tree, "r") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)

    # Load frequencies
    with open(args.frequencies, "r") as fh:
        json_frequencies = json.load(fh)

    kde_frequencies = KdeFrequencies.from_json(json_frequencies)
    start_date = kde_frequencies.start_date
    end_date = kde_frequencies.end_date

    # Setup a model to test LBI

    # The initial model can be configured once and executed several times with
    # different parameters to avoid recalculating censored frequencies, etc. each
    # time.

    predictor_kwargs = {
        "tau": 0.75,
        "time_window": 0.75
    }
    masks_path = os.path.join(augur_path, "builds", "flu", "metadata", "ha_masks.tsv")

    model = FitnessModel(
        tree,
        kde_frequencies,
        ["lbi"],
        epitope_masks_fname=masks_path,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        predictor_kwargs=predictor_kwargs,
        verbose=0,
        step_size=0.5,
        min_freq=0.1
    )

    model.prep_nodes()
    model.calc_node_frequencies()
    model.select_clades_for_fitting()

    # Run model with a range of LBI parameters
    # Create a range of tau and time window values and then build the Cartesian product of all values.

    tau_values = [0.03125, 0.0625, 0.125, 0.25, 0.5, 0.75, 1.0]
    window_values = [0.1, 0.25, 0.75]
    lbi_parameters = [{"tau": tau, "time_window": window} for tau, window in itertools.product(tau_values, window_values)]

    # Given the prepare fitness model, calculate predictors using the different LBI
    # parameters and learn parameters for those parameters.

    results = []
    for parameter_set in lbi_parameters:
        model.predictor_kwargs.update(parameter_set)
        print("Fitting to " + str(parameter_set))

        model.calc_all_predictors()
        model.standardize_predictors()
        model.learn_parameters()

        correlation_null, correlation_raw, correlation_rel = model.get_correlation()
        correlation = correlation_rel[0]
        mcc = get_matthews_correlation_coefficient_for_data_frame(model.pred_vs_true_df)
        results.append({
            "tau": parameter_set["tau"],
            "time_window": parameter_set["time_window"],
            "correlation": correlation,
            "mcc": mcc,
            "param": model.model_params[0]
        })

        print("Tau: %s, window: %s, correlation: %s, MCC: %s, param: %s" % (
                parameter_set["tau"],
                parameter_set["time_window"],
                correlation,
                mcc,
                model.model_params[0]
        ))

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.results, sep="\t", index=False)
