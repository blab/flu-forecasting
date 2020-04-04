import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--titers", required=True, help="JSON of complete titer records to be filtered by passage type")
    parser.add_argument("--passage-type", required=True, help="type of passage for viruses used in titer assays")
    parser.add_argument("--output", required=True, help="table of filtered titer records by passage type")

    args = parser.parse_args()

    df = pd.read_json(args.titers)
    passaged = (df["serum_passage_category"] == args.passage_type)
    tdb_passaged = df["index"].apply(lambda index: isinstance(index, list) and args.passage_type in index)
    tsv_fields = [
        "virus_strain",
        "serum_strain",
        "serum_id",
        "source",
        "titer",
        "assay_type"
    ]

    titers_df = df.loc[(passaged | tdb_passaged), tsv_fields]
    titers_df.to_csv(args.output, sep="\t", header=False, index=False)
