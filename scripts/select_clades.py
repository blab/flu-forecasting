"""Select clades per timepoint to use when training/validating a fitness model.
"""
import argparse
import numpy as np
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Standardize predictors",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tip-attributes", required=True, help="tab-delimited file describing tip attributes at all timepoints")
    parser.add_argument("--tips-to-clades", required=True, help="tab-delimited file of all clades per tip and timepoint")
    parser.add_argument("--delta-months", required=True, type=int, help="number of months to project clade frequencies into the future")
    parser.add_argument("--output", required=True, help="tab-delimited file of clades per timepoint and their corresponding tips and tip frequencies at the given delta time in the future")
    args = parser.parse_args()

    delta_time_offset = pd.DateOffset(months=args.delta_months)

    # Load tip attributes, subsetting to relevant frequency and time information.
    tips = pd.read_csv(args.tip_attributes, sep="\t", parse_dates=["timepoint"])
    tips = tips.loc[:, ["strain", "clade_membership", "timepoint", "frequency"]].copy()

    # Confirm tip frequencies sum to 1 per timepoint.
    summed_tip_frequencies = tips.groupby("timepoint")["frequency"].sum()
    print(summed_tip_frequencies)
    assert all([
        np.isclose(total, 1.0, atol=1e-3)
        for total in summed_tip_frequencies
    ])

    # Identify distinct clades per timepoint.
    clades = tips.loc[:, ["timepoint", "clade_membership"]].drop_duplicates().copy()
    clades = clades.rename(columns={"timepoint": "initial_timepoint"})

    # Annotate future timepoint.
    clades["final_timepoint"] = clades["initial_timepoint"] + delta_time_offset

    # Load mapping of tips to all possible clades at each timepoint.
    tips_to_clades = pd.read_csv(args.tips_to_clades, sep="\t", parse_dates=["timepoint"])
    tips_to_clades = tips_to_clades.loc[:, ["tip", "clade_membership", "depth", "timepoint"]].copy()

    # Get all tip-clade combinations by timepoint for the distinct clades.
    future_tips_by_clades = clades.merge(
        tips_to_clades,
        how="inner",
        left_on=["final_timepoint", "clade_membership"],
        right_on=["timepoint", "clade_membership"]
    )

    # Drop redundant columns.
    future_tips_by_clades = future_tips_by_clades.drop(
        columns=["timepoint"]
    )

    # Get the closest clade to each tip by timepoint. This relies on records
    # being sorted by depth of clade from tip.
    future_tips_by_clades = future_tips_by_clades.sort_values(["initial_timepoint", "tip", "depth"]).groupby(["initial_timepoint", "tip"]).first().reset_index()

    # Get frequencies of future tips associated with current clades.
    future_clade_frequencies = future_tips_by_clades.merge(tips, how="inner", left_on=["tip", "final_timepoint"], right_on=["strain", "timepoint"], suffixes=["", "_tip"])
    future_clade_frequencies = future_clade_frequencies.drop(
        columns=[
            "tip",
            "depth",
            "clade_membership_tip",
            "timepoint"
        ]
    )

    # Confirm that future frequencies sum to 1.
    print(future_clade_frequencies.groupby("initial_timepoint")["frequency"].sum())

    # Confirm the future frequencies of individual clades.
    print(future_clade_frequencies.groupby(["initial_timepoint", "clade_membership"])["frequency"].sum())

    # Left join original clades table with the future tip frequencies to enable
    # assessment of all current clades including those without future tips.
    final_clade_frequencies = clades.merge(
        future_clade_frequencies,
        how="left",
        on=["initial_timepoint", "final_timepoint", "clade_membership"]
    )

    # Fill frequency of clades without any future tips with zeros to enable a
    # simple groupby in the future to get observed future frequencies of all
    # clades.
    final_clade_frequencies["frequency"] = final_clade_frequencies["frequency"].fillna(0.0)

    # Confirm that future frequencies sum to 1.
    print(final_clade_frequencies.groupby("initial_timepoint")["frequency"].sum())

    # Save clade future tip frequencies by timepoint.
    final_clade_frequencies.to_csv(args.output, sep="\t", na_rep="N/A", index=False)
