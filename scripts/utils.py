"""
Utilities for working with augur.
"""
import json
import os
import sys

# augur imports.
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "augur"))
from base.frequencies import KdeFrequencies
from base.io_util import json_to_tree


def load_tree_from_json_filename(filename):
    # Load JSON tree.
    with open(filename, "r") as json_fh:
         json_tree = json.load(json_fh)

    # Convert JSON tree layout to a Biopython Clade instance.
    tree = json_to_tree(json_tree)

    return tree


def load_frequencies_from_json_filename(filename):
    # Load JSON frequencies.
    with open(filename, "r") as json_fh:
         json_frequencies = json.load(json_fh)

    # Create frequencies instance from JSON.
    frequencies = KdeFrequencies.from_json(json_frequencies)

    return frequencies
