import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tip-attributes", required=True, help="table of tip attributes from one or more timepoints")
    parser.add_argument("--output", required=True, help="table of tip attributes annotated with a 'naive' predictor")

    args = parser.parse_args()

    # Annotate a predictor for a naive model with no growth.
    df = pd.read_csv(args.tip_attributes, sep="\t")
    df["naive"] = 0.0
    df.to_csv(args.output, sep="\t", index=False)
