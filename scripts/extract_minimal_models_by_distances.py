import argparse
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--model", required=True, help="JSON for a complete fitness model output")
    parser.add_argument("--output", required=True, help="JSON for a minimal fitness model (coefficients only)")

    args = parser.parse_args()

    with open(args.model, "r") as fh:
        model = json.load(fh)

    if "scores" in model:
        del model["scores"]

    with open(args.output, "w") as oh:
        json.dump(model, oh, indent=1)
