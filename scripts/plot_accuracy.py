import argparse
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Annotate predictor class
def annotate_predictor_class(predictor):
    if predictor in ["ep", "ep-ep_x", "na_ep", "cTiter", "cTiterSub"]:
        pclass = "antigenic"
    elif predictor in ["ne", "dms"]:
        pclass = "functional"
    elif predictor in ["lbi"]:
        pclass = "phylogenetic"
    elif predictor == "all" or predictor in ["ep-ep_x-ne", "ep-ep_x-ne-na_ep", "cTiterSub-lbi", "ne-cTiterSub-lbi", "ep-ep_x-ne-cTiterSub-lbi"]:
        pclass = "ensemble"
    else:
        pclass = "control"

    return pclass

# Annotate predictor order
def annotate_order(predictor):
    predictor_order = {
        "null": 0,
        "ep": 1,
        "ep-ep_x": 2,
        "na_ep": 2.5,
        "cTiterSub": 3,
        "cTiter": 3.5,
        "ne": 4,
        "dms": 4,
        "lbi": 5,
        "ep-ep_x-ne": 6,
        "ep-ep_x-ne-na_ep": 6.4,
        "cTiterSub-lbi": 6.5,
        "ne-cTiterSub-lbi": 6.6,
        "ep-ep_x-ne-cTiterSub-lbi": 6.7,
        "all": 7
    }

    return predictor_order.get(predictor, 100)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("accuracy")
    parser.add_argument("correlation")
    parser.add_argument("mcc")

    args = parser.parse_args()

    fontsize = 18
    params = {
        'axes.labelsize': fontsize,
        'font.size': fontsize,
        'legend.fontsize': 16,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize,
        'text.usetex': False,
        'figure.figsize': [8, 6],
        'savefig.dpi': 300,
        'figure.dpi': 120
    }
    plt.rcParams.update(params)

    accuracy_df = pd.read_table(args.accuracy, keep_default_na=False, na_values=["NaN"])
    accuracy_df["predictor_class"] = accuracy_df["predictors"].apply(annotate_predictor_class)
    accuracy_df["order"] = accuracy_df["predictors"].apply(annotate_order)
    accuracy_df["build"] = accuracy_df.apply(lambda row: "%s, %sv" % (row["year_range"], row["viruses"]), axis=1)

    accuracy_df = accuracy_df.sort_values(["order"])

    g = sns.catplot(
        y="predictors",
        x="correlation_rel",
        row="build",
        data=accuracy_df,
        kind="point",
        ci="sd",
        join=False,
        height=10,
        aspect=1.33,
        orient="h"
    )
    g.set_axis_labels("Frequency correlation", "Predictors")
    plt.savefig(args.correlation)

    g = sns.catplot(
        y="predictors",
        x="mcc",
        row="build",
        data=accuracy_df,
        kind="point",
        ci="sd",
        join=False,
        height=10,
        aspect=1.33,
        orient="h"
    )
    g.set_axis_labels("Matthew's correlation coefficient", "Predictors")
    plt.savefig(args.mcc)
