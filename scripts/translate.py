"""Translate nucleotide sequences to amino acid sequences for the requested genes.
"""
import argparse
from augur.translate import translate_feature
from augur.utils import load_features
import Bio.Seq
import Bio.SeqIO


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Translate nucleotide sequences to amino acid sequences for the requested genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--sequences', required=True, help='FASTA file of nucleotide sequences to translate')
    parser.add_argument('--reference-sequence', required=True, help='GenBank or GFF file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to translate (list or file containing list)")
    parser.add_argument('--output', nargs='+', help="FASTA files of amino acid sequences per gene")

    args = parser.parse_args()

    # Load features for requested genes.
    features = load_features(args.reference_sequence, args.genes)

    # Load sequences indexed by sequence id.
    sequences = {
        sequence.id: str(sequence.seq)
        for sequence in Bio.SeqIO.parse(args.sequences, "fasta")
        if "N" not in str(sequence.seq)
    }
    #if sorted(set(list(str(sequence.seq)))) == ["A", "C", "G", "T"]

    # Translate requested genes.
    translations = {}
    invalid_samples = set()
    for feature_name, output_file in zip(args.genes, args.output):
        translations[feature_name] = translate_feature(sequences, features[feature_name])
        records = [
            Bio.SeqIO.SeqRecord(
                Bio.Seq.Seq(sequence),
                id=sample,
                name="",
                description=""
            )
            for sample, sequence in translations[feature_name].items()
            if sample not in invalid_samples
        ]
        print("Writing %i records to %s" % (len(records), output_file))
        Bio.SeqIO.write(records, output_file, "fasta")
