"""Convert a given titer substitution model JSON to a distance map supported by augur distance.
"""
import argparse
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert titer substitution model to distance map",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--model", required=True, help="JSON from titer substitution model")
    parser.add_argument("--output", required=True, help="distance map JSON")
    args = parser.parse_args()

    # Load titer model.
    with open(args.model, "r") as fh:
        model = json.load(fh)

    # Prepare a distance map for the model.
    distance_map = {
        "name": "titer_substitution_model",
        "default": 0.0,
        "map": {}
    }

    # Convert values like:
    # "HA1:E173K": 0.4656
    # to distance map format.
    for substitution, weight in model["substitution"].items():
        gene, mutation = substitution.split(":")
        ancestral = mutation[0]
        derived = mutation[-1]
        position = mutation[1:-1]

        if ancestral != "X" and derived != "X":
            if gene not in distance_map["map"]:
                distance_map["map"][gene] = {}

            if position not in distance_map["map"][gene]:
                distance_map["map"][gene][position] = []

            distance_map["map"][gene][position].append({
                "from": ancestral,
                "to": derived,
                "weight": weight
            })

    # Save the distance map.
    with open(args.output, "w") as oh:
        json.dump(distance_map, oh, sort_keys=True, indent=1)
