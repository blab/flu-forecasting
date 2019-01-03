"""Plot fitness on a phylogenetic tree.
"""
import argparse
import Bio
import Bio.Phylo
import json
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.collections import LineCollection, PatchCollection
import matplotlib.patches as mpatches
import numpy as np
import os
import sys

# augur imports.
augur_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "augur")
sys.path.append(augur_path)
from base.fitness_model import fitness_model as FitnessModel
from base.frequencies import KdeFrequencies
from base.io_util import json_to_tree


def plot_tree(tree, model, timepoint, color_by_trait, initial_branch_width, tip_size, figure_name=None, min_date=None, limit_max_date=False):
    """Plot a BioPython Phylo tree in the BALTIC-style.
    """
    # Tips being fit at this timepoint by clade id.
    tips_being_fit = set()
    for clade in model.fit_clades[timepoint]:
        for tip in clade.get_terminals():
            tips_being_fit.add(tip.clade)

    # Plot H3N2 tree in BALTIC style from Bio.Phylo tree.
    mpl.rcParams['savefig.dpi'] = 120
    mpl.rcParams['figure.dpi'] = 100

    mpl.rcParams['font.weight']=300
    mpl.rcParams['axes.labelweight']=300
    mpl.rcParams['font.size']=14

    yvalues = [node.yvalue for node in tree.find_clades()]
    y_span = max(yvalues)
    y_unit = y_span / float(len(yvalues))

    # Setup colors.
    trait_name = color_by_trait
    traits = [k.attr[trait_name] for k in tree.get_terminals()]
    norm = mpl.colors.Normalize(-2 * np.std(traits), 2 * np.std(traits))
    cmap = mpl.cm.viridis

    #
    # Setup the figure grid.
    #

    fig = plt.figure(figsize=(30, 22), facecolor='w')
    gs = gridspec.GridSpec(2, 1, height_ratios=[30, 1], width_ratios=[1], hspace=0.1, wspace=0.1)
    ax = fig.add_subplot(gs[0])
    colorbar_ax = fig.add_subplot(gs[1])

    L=len([k for k in tree.find_clades() if k.is_terminal()])

    # Setup arrays for tip and internal node coordinates.
    tip_circles_x = []
    tip_circles_y = []
    tip_circles_color = []
    tip_circle_sizes = []
    node_circles_x = []
    node_circles_y = []
    node_circles_color = []
    node_line_widths = []
    node_line_segments = []
    node_line_colors = []
    branch_line_segments = []
    branch_line_widths = []
    branch_line_colors = []
    branch_line_labels = []

    for k in tree.find_clades(): ## iterate over objects in tree
        x=k.attr["num_date"] ## or from x position determined earlier
        y=y_span - k.yvalue ## get y position from .drawTree that was run earlier, but could be anything else

        if k.parent is None:
            xp = None
        else:
            xp=k.parent.attr["num_date"] ## get x position of current object's parent

        if x==None: ## matplotlib won't plot Nones, like root
            x=0.0
        if xp==None:
            xp=x

        c = (0.7, 0.7, 0.7)
        # Only show colors for nodes sampled prior to the current timepoint.
        # Only show colors for nodes if they were in a clade that was being fit
        # at this timepoint.
        if k.attr["num_date"] <= timepoint and k.attr.has_key(trait_name) and k.clade in tips_being_fit:
            c = cmap(norm(k.attr[trait_name]))

        branchWidth=2
        if k.is_terminal(): ## if leaf...
            s = tip_size ## tip size can be fixed

            tip_circle_sizes.append(s)
            tip_circles_x.append(x)
            tip_circles_y.append(y)
            tip_circles_color.append(c)
        else: ## if node...
            k_leaves = [child
                        for child in k.find_clades()
                        if child.is_terminal()]

            # Scale branch widths by the number of tips.
            branchWidth += initial_branch_width * len(k_leaves) / float(L)

            if len(k.clades)==1:
                node_circles_x.append(x)
                node_circles_y.append(y)
                node_circles_color.append(c)

            ax.plot([x,x],[y_span - k.clades[-1].yvalue, y_span - k.clades[0].yvalue], lw=branchWidth, color="#cccccc", ls='-', zorder=9, solid_capstyle='round')

        branch_line_segments.append([(xp, y), (x, y)])
        branch_line_widths.append(branchWidth)

        # Highlight branches of clades that are being fit at the given timepoint.
        if k in model.fit_clades[timepoint]:
            branch_line_colors.append("#00ff00")
        else:
            branch_line_colors.append("#cccccc")

    branch_lc = LineCollection(branch_line_segments, zorder=9)
    branch_lc.set_color(branch_line_colors)
    branch_lc.set_linewidth(branch_line_widths)
    branch_lc.set_label(branch_line_labels)
    branch_lc.set_linestyle("-")
    ax.add_collection(branch_lc)

    # Add circles for tips and internal nodes.
    tip_circle_sizes = np.array(tip_circle_sizes)
    ax.scatter(tip_circles_x, tip_circles_y, s=tip_circle_sizes, facecolor=tip_circles_color, edgecolor='none',zorder=11) ## plot circle for every tip
    ax.scatter(tip_circles_x, tip_circles_y, s=tip_circle_sizes*2, facecolor='k', edgecolor='none', zorder=10) ## plot black circle underneath
    ax.scatter(node_circles_x, node_circles_y, facecolor=node_circles_color, s=50, edgecolor='none', zorder=10, lw=2, marker='|') ## mark every node in the tree to highlight that it's a multitype tree

    # Grey out the time after the current timepoint + delta.
    delta_timepoint = timepoint + model.delta_time
    rectangle = plt.Rectangle((delta_timepoint, 0), model.timepoints[-1] - delta_timepoint, y_span, fc='#999999',
                              alpha=0.2)
    ax.add_patch(rectangle)

    # Highlight the time between the current timepoint and the future timepoint.
    rectangle = plt.Rectangle((timepoint, 0), model.delta_time, y_span, fc='#ff0000',
                              alpha=0.2)
    ax.add_patch(rectangle)

    ax.spines['top'].set_visible(False) ## no axes
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.set_xticks(model.timepoints)
    ax.grid(axis='x',ls='-',color='grey')
    ax.tick_params(axis='y',size=0)
    ax.set_yticklabels([])

    if min_date:
        ax.set_xlim(left=min_date)

    if limit_max_date:
        ax.set_xlim(right=delta_timepoint)
        max_y = max([y_span - node.yvalue for node in tree.find_clades() if node.attr["num_date"] <= delta_timepoint])
        ax.set_ylim(top=max_y)

    ax.set_title("Timepoint: %s to %s" % (timepoint, delta_timepoint))

    cb1 = mpl.colorbar.ColorbarBase(
        colorbar_ax,
        cmap=cmap,
        norm=norm,
        orientation='horizontal'
    )
    cb1.set_label(color_by_trait)

    gs.tight_layout(fig)

    if figure_name:
        plt.savefig(figure_name)
        plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot fitness per tip per timepoint.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("tree", help="JSON for an auspice tree with annotations for fitness model")
    parser.add_argument("frequencies", help="JSON from frequencies estimated for the given tree")
    parser.add_argument("model", help="JSON for fitness model")
    parser.add_argument("output_prefix", help="prefix for output filenames")
    parser.add_argument("--min-date", type=float, help="min date to display on x-axis of tree")
    parser.add_argument("--limit-max-date", action="store_true", help="limit max date to display on x-axis to the end of the projected time")

    args = parser.parse_args()

    # Load tree.
    with open(args.tree, "r") as fh:
        json_tree = json.load(fh)

    tree = json_to_tree(json_tree)

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        json_frequencies = json.load(fh)

    kde_frequencies = KdeFrequencies.from_json(json_frequencies)

    # Load fitness model.
    with open(args.model, "r") as fh:
        json_model = json.load(fh)

    predictors = {record["predictor"]: [round(record["param"], 2), round(record["global_sd"], 2)]
                  for record in json_model["params"]}
    predictors_key = "-".join(sorted([record["predictor"] for record in json_model["params"]]))
    predictor_kwargs = json_model["predictor_kwargs"]

    model = FitnessModel(
        tree,
        kde_frequencies,
        predictors,
        epitope_masks_fname="%s/builds/flu/metadata/ha_masks.tsv" % augur_path,
        epitope_mask_version="wolf",
        tolerance_mask_version="HA1",
        min_freq=0.1,
        predictor_kwargs=predictor_kwargs,
        step_size=json_model["step_size"]
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

    model.select_nonoverlapping_clades_for_fitting()

    timepoints = model.timepoints[:-1]
    for timepoint in timepoints:
        print("Plotting %s" % timepoint)
        # Assign standardized predictor fitness value.
        for tip in tree.get_terminals():
            index = model.tips.index(tip)
            tip.attr["fitness"] = model.fitness(model.model_params, model.predictor_arrays[timepoint][index])

        plot_tree(
            tree,
            model,
            timepoint,
            color_by_trait="fitness",
            initial_branch_width=30.0,
            tip_size=100.0,
            figure_name="%s_%s_%s.png" % (args.output_prefix, "-".join(model.predictors), timepoint),
            min_date=args.min_date,
            limit_max_date=args.limit_max_date
        )
