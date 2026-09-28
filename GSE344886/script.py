# I want to run the same script as R in the python environment

import numpy as numpy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

# Loading counts

counts = pd.read_csv("./GSE344886_GEO_raw_count_6samples.tsv", sep = "\t", index_col = 0)

sample_map = {
    "BEV_1": "GSM9986991",
    "BEV_2": "GSM9986992",
    "BEV_3": "GSM9986993",
    "CTRL_1": "GSM9986994",
    "CTRL_2": "GSM9986995",
    "CTRL_3": "GSM9986996"
}

counts = counts.rename(columns=sample_map) # with the renamed columns having GSE ids

# Loading metadata

all_metadata = {}

with open("GSE344886_series_matrix.txt", "r") as f:
    for line in f:
        if line.startswith("!Sample_") and "\t" in line:
            key, value = line.rstrip("\n").split("\t", 1)
            all_metadata[key] = value.split("\t")

all_metadata = pd.DataFrame(all_metadata)

# Editing the metadata column only for the columns
metadata = all_metadata[["!Sample_geo_accession", "!Sample_characteristics_ch1"]].copy()

rename_mapping = {
    "!Sample_geo_accession": "geo",
    "!Sample_characteristics_ch1": "condition"
}
metadata = metadata.rename(columns=rename_mapping)

metadata["geo"] = metadata["geo"].str.replace('"', '')
metadata["condition"] = metadata["condition"].str.replace("treatment: ", "").str.strip()

metadata["geo"] = metadata["geo"].astype("category")
metadata["condition"] = metadata["condition"].astype("category") # Converting metadata category into category

metadata = metadata.set_index("geo") # Changing Index