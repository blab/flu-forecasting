#!/bin/bash

java -jar ~/projects/santa-sim/dist/santa.jar influenza_h3n2_ha_epochs.xml

seqtk sample simulated_HA_sequences.fasta 0.05 > simulated_HA_sequences_5pct.fasta
#seqtk sample simulated_HA_sequences.fasta 0.2 > simulated_HA_sequences_20pct.fasta

augur parse --sequences simulated_HA_sequences_5pct.fasta --output-sequences sequences.fasta --output-metadata original_metadata.tsv --fields strain generation fitness
#augur parse --sequences simulated_HA_sequences_20pct.fasta --output-sequences sequences.fasta --output-metadata original_metadata.tsv --fields strain generation

# Manual python step to add dates to metadata from generations
# from augur.frequency_estimators import float_to_datestring
# import pandas as pd
# df = pd.read_table("original_metadata.tsv")
# df["num_date"] = 2000.0 + (df["generation"] / 100.0)
# df["date"] = df["num_date"].apply(float_to_datestring)
# df.to_csv("metadata.tsv", sep="\t", header=True, index=False)
python3 correct_metadata.py

# Skip alignment for simulated sequences
augur tree --alignment sequences.fasta --output tree_raw.nwk --nthreads 4

augur refine --alignment sequences.fasta --tree tree_raw.nwk --metadata metadata.tsv --output-tree tree.nwk --output-node-data branch_lengths.json --timetree --clock-filter-iqd 4 --date-inference marginal --date-confidence --coalescent const

augur ancestral --tree tree.nwk --alignment sequences.fasta --inference joint --output nt_muts.json

augur translate --tree tree.nwk --ancestral-sequences nt_muts.json --reference-sequence ../../config/reference_h3n2_ha.gb --output aa_muts.json

augur lbi --tree tree.nwk --branch-lengths branch_lengths.json --output lbi.json --attribute-names lbi --tau 0.25 --window 0.75

augur reconstruct-sequences --tree tree.nwk --gene SigPep --mutations aa_muts.json --internal-nodes --output aa-seq_SigPep.fasta
augur reconstruct-sequences --tree tree.nwk --gene HA1 --mutations aa_muts.json --internal-nodes --output aa-seq_HA1.fasta
augur reconstruct-sequences --tree tree.nwk --gene HA2 --mutations aa_muts.json --internal-nodes --output aa-seq_HA2.fasta

augur distance --tree tree.nwk --alignment aa-seq_SigPep.fasta aa-seq_HA1.fasta aa-seq_HA2.fasta --gene-names SigPep HA1 HA2 --attribute-name ep ne --compare-to root root --map ../../config/distance_maps/h3n2/ha/wolf.json ../../config/distance_maps/h3n2/ha/wolf_nonepitope.json --output distances.json

augur export --tree tree.nwk --metadata metadata.tsv --auspice-config ../../config/auspice_config.json --colors ../../config/colors.tsv --output-tree test_tree.json --output-meta test_meta.json --node-data branch_lengths.json nt_muts.json aa_muts.json lbi.json distances.json --output-sequence test_seq.json

augur frequencies --method kde --metadata metadata.tsv --tree tree.nwk --output test_tip-frequencies.json
