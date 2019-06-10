"""Calculate weighted distances between samples in a given timepoint and both other samples in that timepoint and samples from a timepoint at a given delta time in the future.
"""
import argparse
import numpy as np
import pandas as pd
import sys


def get_distances_by_sample_names(distances):
    """Return a dictionary of distances by pairs of sample names.

    Parameters
    ----------
    distances : pandas.DataFrame
        data frame with the columns distance, sample, and other_sample

    Returns
    -------
    dict :
        dictionary of distances by pairs of sample names
    """
    distances_by_sample_names = {}
    for distance, sample_b, sample_a in distances.values:
        if sample_a not in distances_by_sample_names:
            distances_by_sample_names[sample_a] = {}

        if sample_b not in distances_by_sample_names:
            distances_by_sample_names[sample_b] = {}

        distances_by_sample_names[sample_a][sample_b] = distance
        distances_by_sample_names[sample_b][sample_a] = distance

    return distances_by_sample_names


def get_distance_matrix_by_sample_names(samples_a, samples_b, distances):
    """Return a matrix of distances between pairs of given sample sets.

    Parameters
    ----------
    samples_a, samples_b : list
        names of samples whose pairwise distances should populate the matrix
        with the first samples in rows and the second samples in columns

    distances : dict
        dictionary of distances by pairs of sample names

    Returns
    -------
    ndarray :
        matrix of pairwise distances between the given samples


    >>> samples_a = ["a", "b"]
    >>> samples_b = ["c", "d"]
    >>> distances = {"a": {"c": 1, "d": 2}, "b": {"c": 3, "d": 4}}
    >>> get_distance_matrix_by_sample_names(samples_a, samples_b, distances)
    array([[1., 2.],
           [3., 4.]])
    >>>

    """
    matrix = np.zeros((len(samples_a), len(samples_b)))
    for i, sample_a in enumerate(samples_a):
        for j, sample_b in enumerate(samples_b):
            matrix[i, j] = distances[sample_a][sample_b]

    return matrix


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
    distances_by_sample_names = get_distances_by_sample_names(distances)

    # Find valid timepoints for calculating distances to the future.
    timepoints = tips["timepoint"].drop_duplicates()
    last_timepoint = timepoints.max() - pd.DateOffset(months=args.delta_months)
    valid_timepoints = timepoints[timepoints <= last_timepoint]

    # Calculate weighted distance to the present and future for each sample at a
    # given timepoint.
    weighted_distances = []
    for timepoint in valid_timepoints:
        future_timepoint = timepoint + pd.DateOffset(months=args.delta_months)
        timepoint_tips = tips[tips["timepoint"] == timepoint]
        future_timepoint_tips = tips[tips["timepoint"] == future_timepoint]

        for current_tip, current_tip_frequency in timepoint_tips.loc[:, ["strain", "frequency"]].values:
            weighted_distance_to_present = 0.0
            for other_current_tip, other_current_tip_frequency in timepoint_tips.loc[:, ["strain", "frequency"]].values:
                weighted_distance_to_present += other_current_tip_frequency * distances_by_sample_names[current_tip][other_current_tip]

            weighted_distance_to_future = 0.0
            for future_tip, future_tip_frequency in future_timepoint_tips.loc[:, ["strain", "frequency"]].values:
                weighted_distance_to_future += future_tip_frequency * distances_by_sample_names[current_tip][future_tip]

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
