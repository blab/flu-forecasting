import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata", required=True, help="table of metadata to be filtered based on a date column")
    parser.add_argument("--strains", required=True, help="text file with one strain per line that should be included in the output")
    parser.add_argument("--output", required=True, help="table of filtered metadata")

    args = parser.parse_args()

    metadata = pd.read_table(args.metadata)
    strains = pd.read_table(args.strains, header=None, names=["strain"])

    selected_metadata = strains.merge(metadata, how="left", on="strain")
    selected_metadata.to_csv(args.output, sep="\t", index=False)
