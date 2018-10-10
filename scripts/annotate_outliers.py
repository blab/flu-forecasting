import argparse
import pandas as pd
import re


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_results", help="tabular BLAST results of local H3 sequences against NCBI database")
    parser.add_argument("non_h3_strains", help="text file containing a list of non-H3 strains based on NCBI results")
    args = parser.parse_args()

    df = pd.read_table(args.blast_results, names=("qseqid", "sseqid", "pident", "evalue", "bitscore", "stitle"), header=None)
    df["strain"] = df["qseqid"].apply(lambda name: name.split("|")[0])

    df = df.sort_values(["strain", "bitscore"], ascending=[True, False])
    best_df = df.groupby(["strain"]).first()

    # Annotate organism or location of best strain.
    def annotate_location(title):
        match = re.search("A\/([ \w'\(-]+)\/", title)
        if match:
            return match.groups()[0]
        else:
            return "unmatched"

    best_df["location"] = best_df["stitle"].apply(annotate_location)

    # Annotate H3 status.
    best_df["is_h3"] = best_df["stitle"].str.contains("H3")

    # Annotate subtype of search sequence.
    def annotate_subtype(record):
        if not record["is_h3"]:
            return "non-H3"
        elif record["location"] == "swine":
            return "swine"
        else:
            return "H3"

    best_df["subtype"] = best_df.apply(annotate_subtype, axis=1)
    best_df = best_df.reset_index()

    # Export annotations of strains.
    annotations = best_df.loc[:, ["strain", "subtype"]]

    # Summarize annotations.
    print(annotations["subtype"].value_counts())

    # Save a list of strains that are not annotated as H3.
    best_df.loc[best_df["subtype"] != "H3", "strain"].to_csv(
        args.non_h3_strains,
        sep="\t",
        header=False,
        index=False
    )
