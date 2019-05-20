"""
Convert a distances JSON into a single table of values.
"""
import argparse
import json
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert distances JSON to a data frame",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--json", required=True, help="distances JSON")
    parser.add_argument("--output", required=True, help="tab-delimited file with frequency per node at the last available timepoint")
    parser.add_argument("--distance-attribute", required=True, help="name of the distance attribute to export")
    args = parser.parse_args()

    # Load distances.
    with open(args.json, "r") as fh:
        distances_json = json.load(fh)

    distances = distances_json["nodes"]

    # Create one record for each pairwise distance.
    records = []

    for sample, sample_distances in distances.items():
        for other_sample, distance in sample_distances[args.distance_attribute].items():
            records.append({
                "sample": sample,
                "other_sample": other_sample,
                "distance": distance
            })

    # Convert records into a data frame.
    df = pd.DataFrame(records)

    # Save the table.
    df.to_csv(args.output, sep="\t", index=False, header=True)
