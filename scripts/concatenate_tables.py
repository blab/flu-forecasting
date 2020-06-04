"""Concatenate two or more tables as data frames.
"""
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", help="tables to concatenate")
    parser.add_argument("--separator", default="\t", help="separator between columns in the given tables")
    parser.add_argument("--output", help="concatenated table")

    args = parser.parse_args()

    # Concatenate tables.
    df = pd.concat([
        pd.read_csv(table_file, sep=args.separator)
        for table_file in args.tables
    ], ignore_index=True, sort=True)

    df.to_csv(args.output, sep=args.separator, header=True, index=False)
