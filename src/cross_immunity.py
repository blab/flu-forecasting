"""Calculate cross-immunity for individual samples against past samples based on given pairwise distances and scaled by past frequencies.
"""
import argparse
import json

from forecast.fitness_predictors import inverse_cross_immunity_amplitude


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate cross-immunity",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--frequencies", required=True, help="JSON of frequencies per sample")
    parser.add_argument("--distances", required=True, help="JSON of distances between samples")
    parser.add_argument("--distance-attributes", nargs="+", required=True, help="names of attributes to use from the given distances JSON")
    parser.add_argument("--immunity-attributes", nargs="+", required=True, help="names of attributes to use for the calculated cross-immunities")
    parser.add_argument("--decay-factors", nargs="+", required=True, type=float, help="number of months to project clade frequencies into the future")
    parser.add_argument("--output", required=True, help="cross-immunities calculated from the given distances and frequencies")
    args = parser.parse_args()

    # Load frequencies.
    with open(args.frequencies, "r") as fh:
        frequencies = json.load(fh)

    # Identify maximum frequency per sample.
    max_frequency_per_sample = {
        sample: float(max(sample_frequencies))
        for sample, sample_frequencies in frequencies["data"]["frequencies"].items()
    }

    # Load distances.
    with open(args.distances, "r") as fh:
        distances = json.load(fh)

    distances = distances["nodes"]

    """
  "A/Acre/15093/2010": {
   "ep": 9,
   "ne": 8,
   "rb": 3
  },
    """
    # Calculate cross-immunity for distances defined by the given attributes.
    cross_immunities = {}
    for sample, sample_distances in distances.items():
        for distance_attribute, immunity_attribute, decay_factor in zip(args.distance_attributes, args.immunity_attributes, args.decay_factors):
            if distance_attribute not in sample_distances:
                continue

            if sample not in cross_immunities:
                cross_immunities[sample] = {}

            # Calculate inverse cross-immunity amplitude once from all distances
            # to the current sample. This is an increasingly positive value for
            # samples that are increasingly distant from previous samples.
            cross_immunity = 0.0
            for past_sample, distance in sample_distances[distance_attribute].items():
                cross_immunity += max_frequency_per_sample[past_sample] * inverse_cross_immunity_amplitude(
                    distance,
                    decay_factor
                )

            cross_immunities[sample][immunity_attribute] = cross_immunity

    # Export cross-immunities to JSON.
    with open(args.output, "w") as oh:
        json.dump({"nodes": cross_immunities}, oh, indent=1, sort_keys=True)
