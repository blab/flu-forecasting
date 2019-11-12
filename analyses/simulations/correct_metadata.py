from augur.frequency_estimators import float_to_datestring
import pandas as pd

df = pd.read_table("original_metadata.tsv")
df["num_date"] = 2000.0 + (df["generation"] / 150.0)
df["date"] = df["num_date"].apply(float_to_datestring)
df.to_csv("metadata.tsv", sep="\t", header=True, index=False)
