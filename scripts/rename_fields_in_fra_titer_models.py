import argparse
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--titers-model", required=True, help="titer model JSON from augur titers tree")
    parser.add_argument("--output", required=True, help="titer model JSON with renamed fields for FRA data")

    args = parser.parse_args()

    with open(args.titers_model, "r") as fh:
        titers_json = json.load(fh)

    for sample in titers_json["nodes"].keys():
        titers_json["nodes"][sample]["fra_cTiter"] = titers_json["nodes"][sample]["cTiter"]
        titers_json["nodes"][sample]["fra_dTiter"] = titers_json["nodes"][sample]["dTiter"]
        del titers_json["nodes"][sample]["cTiter"]
        del titers_json["nodes"][sample]["dTiter"]

    with open(args.output, "w") as oh:
        json.dump(titers_json, oh, indent=1)
