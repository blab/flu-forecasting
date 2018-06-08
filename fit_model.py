import argparse
import datetime
import json
import numpy as np
import os

# Add augur source to the Python path.
import sys
code_directory = os.path.join(os.getcwd(), "dist", "augur")
sys.path.insert(0, code_directory)

# Load augur modules.
from base.fitness_model import fitness_model as FitnessModel, make_pivots
from base.frequencies import KdeFrequencies
from base.io_util import json_to_tree, json_to_clade_frequencies
from base.titer_model import TiterCollection


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice tree JSON")
    parser.add_argument("frequencies", help="auspice frequencies JSON")
    parser.add_argument("model", help="output model JSON")
    parser.add_argument("predictors", nargs="+", help="one or more predictors to build model for")
    parser.add_argument("--titers", help="tab-delimited file of titer measurements")
    parser.add_argument("--sigma", type=float, default=1 / 12.0, help="Bandwidth for KDE frequencies")

    args = parser.parse_args()
    predictor_kwargs = {}

    # Load JSON tree.
    with open(args.tree, "r") as json_fh:
         json_tree = json.load(json_fh)

    # Convert JSON tree layout to a Biopython Clade instance.
    tree = json_to_tree(json_tree)

    with open(args.frequencies, "r") as json_fh:
        json_frequencies = json.load(json_fh)

    # Get pivots from frequencies JSON.
    pivots = np.array(json_frequencies["pivots"])

    # Use empty frequencies to force the fitness model to recalculate
    # frequencies with the KDE approach.
    frequencies = KdeFrequencies.estimate_frequencies_for_tree(
        tree,
        pivots,
        sigmaNarrow=args.sigma,
        proportionWide=0.0
    )

    # Determine the time interval from the pivots defined in the frequencies JSON.
    start_date = pivot_to_date(min(pivots))
    end_date = pivot_to_date(max(pivots))
    time_interval = (
        end_date,
        start_date
    )

    # If titers were provided, load them for the model to use.
    if args.titers:
        predictor_kwargs = {
            "lam_avi": 2.0,
            "lam_pot": 0.3,
            "lam_drop": 2.0,
            "preferences_file": "%s/builds/flu/metadata/2017-12-07-H3N2-preferences-rescaled.csv" % code_directory
        }
        titers, strains, sources = TiterCollection.load_from_file(args.titers)
        predictor_kwargs["titers"] = titers

    # Run the fitness model.
    model = FitnessModel(
        tree,
        frequencies,
        time_interval,
        args.predictors,
        pivots=np.around(pivots, 2),
        epitope_masks_fname="%s/builds/flu/metadata/ha_masks.tsv" % code_directory,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        min_freq=0.1,
        predictor_kwargs=predictor_kwargs,
        sigma=args.sigma
    )
    model.predict()
    model.validate_prediction()

    # Save resulting model.
    model.to_json(args.model)
