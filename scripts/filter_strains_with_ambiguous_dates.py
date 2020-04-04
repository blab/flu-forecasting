import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="table of metadata to be filtered based on a date column")
    parser.add_argument("--date-field", default="date", help="name of date column in the metadata")
    parser.add_argument("--output", required=True, help="table of filtered metadata")

    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep="\t")

    # Exclude strains with ambiguous collection dates.
    df[~df[args.date_field].str.contains("XX")].to_csv(args.output, sep="\t", header=True, index=False)
