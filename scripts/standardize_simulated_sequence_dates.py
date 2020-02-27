"""Concatenate two or more tables as data frames.
"""
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", nargs="+", help="")
    parser.add_argument("--start-year", default=2000.0, type=float, help="year to start simulated dates from")
    parser.add_argument("--generations-per-year", default=200.0, type=float, help="number of generations to map to a single yeasr")
    parser.add_argument("--output", help="metadata with standardized dates and nonzero fitness records")

    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep="\t")
    df["num_date"] = args.start_year + (df["generation"] / args.generations_per_year)
    df["date"] = df["num_date"].apply(float_to_datestring)
    df["year"]  = pd.to_datetime(df["date"]).dt.year
    df["month"]  = pd.to_datetime(df["date"]).dt.month

    # Omit records with a fitness of zero.
    df[df["fitness"] > 0].to_csv(args.output, header=True, index=False, sep="\t")
