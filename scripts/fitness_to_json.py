"""
Convert fitness data from models to a node data JSON
"""
import argparse
from augur.utils import write_json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert fitness data from models to a node data JSON",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints")
    parser.add_argument("--timepoint", help="current timepoint", required=True)
    parser.add_argument("--attribute-names", nargs="+", help="names of attributes for tips to export to node data JSON", required=True)
    parser.add_argument("--output", help="JSON file with fitness information by node", required=True)
    args = parser.parse_args()

    # Load tip attributes
    tips = pd.read_csv(
        args.tip_attributes,
        sep="\t",
        parse_dates=["timepoint"]
    )

    timepoint = pd.to_datetime(args.timepoint)
    data = tips.loc[tips["timepoint"] == timepoint, ["strain"] + args.attribute_names].to_dict(orient="records")
    fitnesses = {}
    for record in data:
        fitnesses[record["strain"]] = {}
        for attribute in args.attribute_names:
            fitnesses[record["strain"]][attribute] = record[attribute]

    # Write out the node annotations.
    write_json({"nodes": fitnesses}, args.output)
