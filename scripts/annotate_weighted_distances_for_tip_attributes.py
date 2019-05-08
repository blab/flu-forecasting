"""Calculate weighted distances between samples in a given timepoint and both other samples in that timepoint and samples from a timepoint at a given delta time in the future.
"""
import argparse
import numpy as np
import pandas as pd
import sys


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Annotated weighted distances between viruses",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tip-attributes", required=True, help="a tab-delimited file describing tip attributes at one or more timepoints")
    parser.add_argument("--distances", required=True, help="tab-delimited output file with pairwise distances between samples")
    parser.add_argument("--delta-months", required=True, type=int, help="number of months to project clade frequencies into the future")
    parser.add_argument("--output", required=True, help="tab-delimited output file with mean and standard deviation used to standardize each predictor")
    args = parser.parse_args()

    # Load tip attributes.
    tips = pd.read_csv(args.tip_attributes, sep="\t", parse_dates=["timepoint"])

    # Load distances.
    distances = pd.read_csv(args.distances, sep="\t")

    # Map distances by sample names.
    distances_by_strain = {}
    for distance, sample_a, sample_b in distances.values:
        if sample_a not in distances_by_strain:
            distances_by_strain[sample_a] = {}

        if sample_b not in distances_by_strain:
            distances_by_strain[sample_b] = {}

        distances_by_strain[sample_a][sample_b] = distance
        distances_by_strain[sample_b][sample_a] = distance


    # Find valid timepoints for calculating distances to the future.
    timepoints = tips["timepoint"].drop_duplicates()
    last_timepoint = timepoints.max() - pd.DateOffset(months=args.delta_months)
    valid_timepoints = timepoints[timepoints <= last_timepoint]

    # Calculate weighted distance to the present and future for each sample at a
    # given timepoint.
    weighted_distances = []
    for timepoint in valid_timepoints:
        timepoint_tips = tips[tips["timepoint"] == timepoint]
        future_timepoint_tips = tips[tips["timepoint"] == (timepoint + pd.DateOffset(years=1))]

        for current_tip, current_tip_frequency in timepoint_tips.loc[:, ["strain", "frequency"]].values:
            weighted_distance_to_present = 0.0
            for other_current_tip, other_current_tip_frequency in timepoint_tips.loc[:, ["strain", "frequency"]].values:
                weighted_distance_to_present += other_current_tip_frequency * distances_by_strain[current_tip][other_current_tip]

            weighted_distance_to_future = 0.0
            for future_tip, future_tip_frequency in future_timepoint_tips.loc[:, ["strain", "frequency"]].values:
                weighted_distance_to_future += future_tip_frequency * distances_by_strain[current_tip][future_tip]

            weighted_distances.append({
                "timepoint": timepoint,
                "strain": current_tip,
                "weighted_distance_to_present": weighted_distance_to_present,
                "weighted_distance_to_future": weighted_distance_to_future
            })

    weighted_distances = pd.DataFrame(weighted_distances)

    # Calculate the magnitude of the difference between future and present
    # distances for each sample.
    weighted_distances["log2_distance_effect"] = np.log2(
        weighted_distances["weighted_distance_to_future"] /
        weighted_distances["weighted_distance_to_present"]
    )

    # Annotate samples with weighted distances.
    annotated_tips = tips.merge(
        weighted_distances,
        how="left",
        on=["strain", "timepoint"]
    )

    # Save the new data frame.
    annotated_tips.to_csv(args.output, sep="\t", index=False, na_rep="N/A")
