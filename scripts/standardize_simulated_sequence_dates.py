"""Standardize simulated sequence dates to human readable format (YYYY-MM-DD).
"""
import argparse
import datetime
import pandas as pd


def float_to_datestring(time):
    """Convert a floating point date from TreeTime `numeric_date` to a date string
    """
    # Extract the year and remainder from the floating point date.
    year = int(time)
    remainder = time - year

    # Calculate the day of the year (out of 365 + 0.25 for leap years).
    tm_yday = int(remainder * 365.25)
    if tm_yday == 0:
        tm_yday = 1

    # Construct a date object from the year and day of the year.
    date = datetime.datetime.strptime("%s-%s" % (year, tm_yday), "%Y-%j")

    # Build the date string with zero-padded months and days.
    date_string = "%s-%.2i-%.2i" % (date.year, date.month, date.day)

    return date_string


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", help="metadata for simulated sequences")
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
