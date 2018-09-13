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
    elif predictor == "all" or predictor in ["ep-ep_x-ne", "ep-ep_x-ne-na_ep", "cTiterSub-lbi", "ep-ep_x-ne-cTiterSub-lbi"]:
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
        "ep-ep_x-ne-cTiterSub-lbi": 6.6,
        "all": 7
    }

    return predictor_order.get(predictor, 100)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("parameters")
    parser.add_argument("output")

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

    df = pd.read_table(args.parameters, keep_default_na=False)
    df["order"] = df["predictors"].apply(annotate_order)
    df["build"] = df.apply(lambda row: "%s, %sv" % (row["year_range"], row["viruses"]), axis=1)

    df = df.sort_values(["order"])

    g = sns.catplot(
        y="predictors",
        x="param",
        hue="predictor",
        row="build",
        data=df,
        kind="bar",
        height=10,
        aspect=1.33
    )
    g.set_axis_labels("Model parameters", "Predictors")
    plt.savefig(args.output)
