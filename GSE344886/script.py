# I want to run the same script as R in the python environment

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

# Folder check
folder = "python_results"
if folder not in os.listdir():
    os.makedirs(folder)

print(15 * "--", "Running Analysis", 15 * "--")
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
metadata["condition"] = metadata["condition"].str.replace('"', '').str.strip()

metadata["geo"] = metadata["geo"].astype("category")
metadata["condition"] = metadata["condition"].astype("category") # Converting metadata category into category

metadata = metadata.set_index("geo") # Changing Index

# Filtering counts
keep = (counts >= 10).sum(axis=1) >= 3
counts_filtered = counts.loc[keep]

# Creating DESeq2 object
dds = DeseqDataSet(
    counts=counts_filtered.T,
    metadata=metadata,
    design="~condition"
)

# Running DESeq2
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "bevacizumab-treated", "untreated control"])
stat_res.summary()
results = stat_res.results_df.copy()
results = results.sort_values("padj")

# Saving results to python_results folder
results.to_csv(os.path.join(folder, "deseq_results.csv"))

# Creating Plots
# 1. MA Plot
"""
We could just use the stat_res.plot_MA() function to do so, but I got to know about it later.
Just to flex, I wrote the entire code to plot it right down here.
"""
results["significant"] = (results["padj"] < 0.05) & (abs(results["log2FoldChange"]) >= 1) # Getting significant genes

plt.figure(figsize=(8,5))
plt.scatter(
    x=results.loc[~results["significant"], "baseMean"],
    y=results.loc[~results["significant"], "log2FoldChange"],
    color="grey",
    s=10,
    label="Non-significant"
) # for non=significant genes

plt.scatter(
    x=results.loc[results["significant"], "baseMean"],
    y=results.loc[results["significant"], "log2FoldChange"],
    color="orange",
    s=10,
    label="Significant",
    alpha=0.5
) # for significant genes

plt.xscale("log") # x-axis in log scale
plt.axhline(0, linestyle="--", color="k")
plt.legend()
plt.xlabel("Mean expression")
plt.ylabel("log2 fold change")
plt.title("MA Plot", fontweight="bold")
plt.savefig(os.path.join(folder, "MA_plot.png"))
plt.close()

# 2. Creating PCA

norm_counts = dds.layers["normed_counts"]
log_counts = np.log2(norm_counts + 1)

pca = PCA(n_components=2)

pca_result = pca.fit_transform(log_counts)

pca_df = pd.DataFrame({
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
    "names": dds.obs_names,
    "condition": metadata["condition"].values
})

# Plotting PCA

plt.figure(figsize=(8,5))
sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="condition", s=100, palette="Set1")
plt.title("PCA", fontweight="bold")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.savefig(os.path.join(folder, "PCA.png"), bbox_inches="tight")
plt.close()

# Another approach using vst
dds.vst_fit()
vst_counts = dds.vst_transform()
vst_df = pd.DataFrame(vst_counts, index=counts_filtered.T.index, columns=counts_filtered.T.columns)
pca_vst = PCA(n_components=2)
pca_values = pca_vst.fit_transform(vst_df)
pcavst_df = pd.DataFrame(pca_values, columns=["PC1", "PC2"])
pcavst_df["condition"] = metadata.loc[counts.T.index, "condition"].values

sns.scatterplot(data=pcavst_df, x="PC1", y="PC2", hue="condition", s=100, palette="Set1")
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.title("PCA_vst", fontweight="bold")
plt.xlabel(f"PC1 ({pca_vst.explained_variance_ratio_[0] * 100:.2f})")
plt.ylabel(f"PC2 ({pca_vst.explained_variance_ratio_[1] * 100:.2f})")
plt.savefig(os.path.join(folder, "PCA_vst.png"), bbox_inches="tight")
plt.close()

# 3. Volcano Plot
vol_df = results.dropna(subset=["log2FoldChange", "padj"]).copy()
vol_df["-log10(padj)"] = -np.log10(vol_df["padj"].clip(lower=1e-300))

plt.figure(figsize=(8,6))

non_sig = ~vol_df["significant"] # for non-significant genes

plt.scatter(
    vol_df.loc[non_sig, "log2FoldChange"],
    vol_df.loc[non_sig, "-log10(padj)"],
    color="grey",
    s=10,
    alpha=0.5,
    label="Non-significant"
)

upreg = (vol_df["significant"]) & (vol_df["log2FoldChange"] > 0)

plt.scatter(
    vol_df.loc[upreg, "log2FoldChange"],
    vol_df.loc[upreg, "-log10(padj)"],
    color="red",
    s=10,
    alpha=0.8,
    label="Up-regulated"
)

downreg = (vol_df["significant"] & (vol_df["log2FoldChange"] < 0))

plt.scatter(
    vol_df.loc[downreg, "log2FoldChange"],
    vol_df.loc[downreg, "-log10(padj)"],
    color="blue",
    s=10,
    alpha=0.8,
    label="Down-regulated"
)

# Threshold
plt.axhline(-np.log10(0.05), linestyle="--", color="k")
plt.axvline(1, linestyle="--", color="k")
plt.axvline(-1, linestyle="--", color="k")

plt.xlabel("log2 fold change")
plt.ylabel("-log10(adjusted p-value)")
plt.title("Volcano Plot", fontweight="bold")
plt.legend()
plt.savefig(os.path.join(folder, "volcano_plot.png"), bbox_inches="tight")
plt.close()

# Saving upregulated and downregulated gene list to csv
up_df = vol_df[upreg]
down_df = vol_df[downreg]

reg_df = pd.concat([up_df, down_df])

reg_df["case"] = ["Up" if x > 0 else "Down" for x in reg_df["log2FoldChange"]]

reg_df = reg_df[["case", "log2FoldChange", "padj", "-log10(padj)"]]

reg_df.to_csv(os.path.join(folder, "regulated_genes.csv"))