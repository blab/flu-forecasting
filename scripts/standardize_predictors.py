"""Standardize predictor values in a given tip attributes data frame using the mean and standard deviation from a fixed interval of training data and output the standardized attributes data frame and a tab-delimited file of the summary statistics used for standardization of each predictor.
"""
import argparse
import pandas as pd
from sklearn.preprocessing import StandardScaler
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Standardize predictors",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tip-attributes", required=True, help="a tab-delimited file describing tip attributes at one or more timepoints")
    parser.add_argument("--standardized-attributes", required=True, help="tab-delimited output file with standardized values for each given predictor")
    parser.add_argument("--statistics", required=True, help="tab-delimited output file with mean and standard deviation used to standardize each predictor")
    parser.add_argument("--start-date", required=True, help="the earliest timepoint to calculate means and standard deviations from (YYYY-MM-DD format)")
    parser.add_argument("--end-date", required=True, help="the latest timepoint to calculate means and standard deviations from (YYYY-MM-DD format)")
    parser.add_argument("--predictors", nargs="+", help="a list of columns names for predictors whose values should be standardized")
    args = parser.parse_args()

    # Load tip attributes.
    df = pd.read_csv(args.tip_attributes, sep="\t")

    # Confirm presence of all requested predictor columns.
    missing_columns = set(args.predictors) - set(df.columns)
    if len(missing_columns) > 0:
        print("Error: Could not find the following columns in the given attributes table:", file=sys.stderr)
        for column in missing_columns:
            print(f"  - {column}", file=sys.stderr)

        sys.exit(1)

    # Confirm that timepoints are defined.
    if not "timepoint" in df.columns:
        print("Error: The given attributes table is missing a 'timepoint' column", file=sys.stderr)
        sys.exit(1)

    # Convert string timepoints to datetime instances.
    df["timepoint"] = pd.to_datetime(df["timepoint"])

    # Confirm availability of timepoints for calculating summary statistics.
    start_date = pd.to_datetime(args.start_date)
    end_date = pd.to_datetime(args.end_date)
    valid_timepoints = (start_date <= df["timepoint"]) & (df["timepoint"] <= end_date)
    if valid_timepoints.sum() == 0:
        print(f"Error: The requested timepoints ({args.start_date} to {args.end_date}) were not found in the given attributes table", file=sys.stderr)
        sys.exit(1)

    # Fill missing predictor values with the sample mean.
    for predictor in args.predictors:
        df[predictor] = df[predictor].fillna(df[predictor].mean())

    # Extract only the predictor values in the requested time interval.
    predictor_training_df = df.loc[valid_timepoints, args.predictors].copy()

    # Create a standard scaler and fit the requested training interval for the
    # requested predictor columns.
    scaler = StandardScaler()
    scaler.fit(predictor_training_df.values)

    # Transform predictor values across the entire data frame.
    predictor_df = df.loc[:, args.predictors].copy()
    standardized_values = scaler.transform(predictor_df.values)

    # Create a new data frame with the standardized values.
    standardized_df = df.copy()
    standardized_df.update(pd.DataFrame(standardized_values, columns=args.predictors))

    # Save the new data frame.
    standardized_df.to_csv(args.standardized_attributes, sep="\t", index=False)

    # Collect the summary statistics from the scaler.
    stats_df = pd.DataFrame({
        "mean": scaler.mean_,
        "var": scaler.var_,
        "predictor": args.predictors
    })

    # Save summary statistics.
    stats_df.to_csv(args.statistics, sep="\t", index=False)
