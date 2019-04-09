"""Calculate mutation frequencies per tip.
"""
import argparse
from augur.reconstruct_sequences import load_alignments
from Bio.Align import MultipleSeqAlignment
import Bio.AlignIO
from collections import defaultdict
import json
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find clades in a tree by distinct amino acid haplotypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--frequencies", required=True, help="mutation frequencies in auspice JSON format")
    parser.add_argument("--tips", required=True, help="tip attributes with strain name and frequency in TSV format")
    parser.add_argument("--gene-names", required=True, nargs="+", help="gene names corresponding to translations provided")
    parser.add_argument("--translations", required=True, nargs="+", help="FASTA file(s) of amino acid sequences per node")
    parser.add_argument("--max-frequency", type=float, default=0.99, help="maximum frequency any single mutation can have for its site to be considered")
    parser.add_argument("--output", required=True, help="table of mutation frequencies by tip")

    args = parser.parse_args()

    # Load frequencies and convert them into a data frame.
    with open(args.frequencies, "r") as fh:
        mutation_frequencies = json.load(fh)

    mutations_df = pd.DataFrame([
        {"site": mutation[:-1], "residue": mutation[-1], "frequency": frequencies[-1]}
        for mutation, frequencies in mutation_frequencies.items()
        if mutation != "pivots" and not mutation.endswith("counts")
    ])

    # Load translations for tips and index them by gene name and tip name.
    translations = load_alignments(args.translations, args.gene_names)
    translations_by_gene_name = {}
    for gene in translations:
        translations_by_gene_name[gene] = {}
        for seq in translations[gene]:
            translations_by_gene_name[gene][seq.name] = str(seq.seq)

    # Find the maximum frequency of any residue at each site.
    max_frequency_by_site = mutations_df.groupby("site")["frequency"].max().reset_index()

    # Find all sites to consider at the current timepoint based on the maximum
    # mutation frequency. Sites with a lower maximum frequency have at least one
    # other residue at sufficient global frequencies.
    sites_to_track = max_frequency_by_site.loc[max_frequency_by_site["frequency"] < args.max_frequency, "site"].values
    print("Tracking %i sites including: %s" % (len(sites_to_track), sites_to_track))

    # Filter mutations based on sites to track.
    mutations_to_track_df = mutations_df[mutations_df["site"].isin(sites_to_track)].copy()
    print("Tracking %i mutations" % len(mutations_to_track_df))

    # Collect all positions to inspect per gene.
    positions_per_gene = defaultdict(list)
    for site in sites_to_track:
        gene, position = site.split(":")
        positions_per_gene[gene].append(int(position))

    # Find all tips (samples) and residues associated with each tracked site.
    samples_and_sites = []
    for gene, positions in positions_per_gene.items():
        for seq in translations[gene]:
            # Skip internal nodes.
            if seq.name.startswith("NODE"):
                continue

            # Store this record's sequence for each site to track.
            for position in positions:
                samples_and_sites.append({
                    "site": "%s:%s" % (gene, position),
                    "residue": translations_by_gene_name[gene][seq.name][position - 1],
                    "strain": seq.name
                })

    samples_and_sites_df = pd.DataFrame(samples_and_sites)
    print("Found %i tips and residues for all tracked sites" % len(samples_and_sites_df))

    # Connect all mutations with their respective tips by site and residue.
    mutations_per_sample_df = mutations_to_track_df.merge(samples_and_sites_df, on=["site", "residue"])

    # Annotate the original mutation as site and residue.
    mutations_per_sample_df["mutation"] = mutations_per_sample_df["site"] + mutations_per_sample_df["residue"]

    # Load tip (sample) attributes for this timepoint which includes frequency per sample.
    samples_df = pd.read_csv(
        args.tips,
        sep="\t",
        usecols=["strain", "timepoint", "frequency"],
        parse_dates=["timepoint"]
    )
    print("Loaded attributes for %i tips" % (len(samples_df)))

    # Annotate mutations per sample with sample frequencies and timepoint.
    mutations_with_sample_frequencies_df = mutations_per_sample_df.merge(
        samples_df,
        how="inner",
        on=["strain"],
        suffixes=["_mutation", ""]
    ).drop(columns=["frequency_mutation"])

    mutations_with_sample_frequencies_df = mutations_with_sample_frequencies_df[mutations_with_sample_frequencies_df["frequency"] > 0].copy()

    print("Total tip frequencies per site:")
    print(mutations_with_sample_frequencies_df.groupby(["timepoint", "site"])["frequency"].sum())

    print("Total mutation frequencies per site:")
    print(mutations_to_track_df.groupby(["site"])["frequency"].sum())

    # Save final table.
    mutations_with_sample_frequencies_df.to_csv(args.output, sep="\t", index=False)
