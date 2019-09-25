import argparse
from augur.utils import json_to_tree, read_tree
import Bio.Phylo
import json
import matplotlib as mpl
mpl.use("Agg")
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import sys


def timestamp_to_float(time):
    """Convert a pandas timestamp to a floating point date.
    """
    return time.year + ((time.month - 1) / 12.0)


def plot_tree(tree, figure_name, color_by_trait, initial_branch_width, tip_size,
              start_date, end_date, include_color_bar):
    """Plot a BioPython Phylo tree in the BALTIC-style.
    """
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
    traits = [k.attr[trait_name] for k in tree.find_clades()]
    norm = mpl.colors.Normalize(min(traits), max(traits))
    cmap = mpl.cm.viridis

    #
    # Setup the figure grid.
    #

    if include_color_bar:
        fig = plt.figure(figsize=(8, 6), facecolor='w')
        gs = gridspec.GridSpec(2, 1, height_ratios=[14, 1], width_ratios=[1], hspace=0.1, wspace=0.1)
        ax = fig.add_subplot(gs[0])
        colorbar_ax = fig.add_subplot(gs[1])
    else:
        fig = plt.figure(figsize=(8, 4), facecolor='w')
        gs = gridspec.GridSpec(1, 1)
        ax = fig.add_subplot(gs[0])

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
        y=k.yvalue ## get y position from .drawTree that was run earlier, but could be anything else

        if k.parent is None:
            xp = None
        else:
            xp=k.parent.attr["num_date"] ## get x position of current object's parent

        if x==None: ## matplotlib won't plot Nones, like root
            x=0.0
        if xp==None:
            xp=x

        c = 'k'
        if trait_name in k.attr:
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

            ax.plot([x,x],[k.clades[-1].yvalue, k.clades[0].yvalue], lw=branchWidth, color=c, ls='-', zorder=9, solid_capstyle='round')

        branch_line_segments.append([(xp, y), (x, y)])
        branch_line_widths.append(branchWidth)
        branch_line_colors.append(c)

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

    #ax.set_ylim(-10, y_span - 300)

    ax.spines['top'].set_visible(False) ## no axes
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.grid(axis='x',ls='-',color='grey')
    ax.tick_params(axis='y',size=0)
    ax.set_yticklabels([])

    if start_date:
        # Always add a buffer to the left edge of the plot so data up to the
        # given end date can be clearly seen.
        ax.set_xlim(left=timestamp_to_float(pd.to_datetime(start_date)) - 2.0)

    if end_date:
        # Always add a buffer of 3 months to the right edge of the plot so data
        # up to the given end date can be clearly seen.
        ax.set_xlim(right=timestamp_to_float(pd.to_datetime(end_date)) + 0.25)

    if include_color_bar:
        cb1 = mpl.colorbar.ColorbarBase(
            colorbar_ax,
            cmap=cmap,
            norm=norm,
            orientation='horizontal'
        )
        cb1.set_label(color_by_trait)

    gs.tight_layout(fig)
    plt.savefig(figure_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice tree JSON or Newick tree")
    parser.add_argument("output", help="plotted tree figure")
    parser.add_argument("--colorby", help="trait in tree to color by", default="num_date")
    parser.add_argument("--branch_width", help="initial branch width", type=int, default=10)
    parser.add_argument("--tip_size", help="tip size", type=int, default=10)
    parser.add_argument("--start-date", help="earliest date to show on the x-axis")
    parser.add_argument("--end-date", help="latest date to show on the x-axis")
    parser.add_argument("--include-color-bar", action="store_true", help="display a color bar for the color by option at the bottom of the plot")
    args = parser.parse_args()

    if args.tree.endswith(".json"):
        with open(args.tree, "r") as json_fh:
            json_tree = json.load(json_fh)

        # Convert JSON tree layout to a Biopython Clade instance.
        tree = json_to_tree(json_tree)

        # Plot the tree.
        plot_tree(
            tree,
            args.output,
            args.colorby,
            args.branch_width,
            args.tip_size,
            args.start_date,
            args.end_date,
            args.include_color_bar
        )
    else:
        tree = read_tree(args.tree)
        tree.ladderize()
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        Bio.Phylo.draw(tree, axes=ax, label_func=lambda node: "", show_confidence=False)
        plt.savefig(args.output)
