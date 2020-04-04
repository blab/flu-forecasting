import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--node-data", required=True, help="table of node data from one or more timepoints")
    parser.add_argument("--kde-frequencies", required=True, help="table of KDE frequencies by strain name and timepoint")
    parser.add_argument("--diffusion-frequencies", required=True, help="table of diffusion frequencies by strain name and timepoint")
    parser.add_argument("--preferred-frequency-method", choices=["kde", "diffusion"], help="specify which frequency method should be used for the primary frequency column")
    parser.add_argument("--output", required=True, help="table of merged node data and frequencies")

    args = parser.parse_args()

    node_data = pd.read_table(args.node_data)
    kde_frequencies = pd.read_table(args.kde_frequencies)
    diffusion_frequencies = pd.read_table(args.diffusion_frequencies)
    df = node_data.merge(
        kde_frequencies,
        how="inner",
        on=["strain", "timepoint", "is_terminal"]
    ).merge(
        diffusion_frequencies,
        how="inner",
        on=["strain", "timepoint", "is_terminal"]
    )

    # Annotate frequency by the preferred method if there isn't already a
    # frequency column defined.
    if "frequency" not in df.columns:
        df["frequency"] = df["%s_frequency" % args.preferred_frequency_method]

    df = df[df["frequency"] > 0.0].copy()
    df.to_csv(args.output, sep="\t", index=False, header=True)
