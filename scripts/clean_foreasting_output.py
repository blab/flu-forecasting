import os
import sys

sample = sys.argv[1]

if "natural" in sample:
    type = "natural"
elif "simulated" in sample:
    type = "simulated"

os.system(f'find ./results/builds/{type}/{sample} -name "node_data.tsv" -delete')
os.system(f'find ./results/builds/{type}/{sample} -name "tip_attributes*.tsv" -delete')
os.system(f'find ./results -name "distance_model*.tsv" -delete')

## Example command: "python clean_run.py natural_sample_1_with_90_vpm_sliding"
