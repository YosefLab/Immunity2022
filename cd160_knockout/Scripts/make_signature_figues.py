#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import feather
import json
import os

# %% Define some colors and styles

COOL_RED = sns.color_palette("coolwarm")[5]
COOL_BLUE = sns.color_palette("coolwarm")[0]
colors = plt.get_cmap("Dark2").colors

sns.set_style("whitegrid")
sns.set_context("notebook")
plt.rcParams["figure.dpi"] = 100
plt.rcParams["ytick.left"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["legend.framealpha"] = 1.0
plt.rcParams["patch.edgecolor"] = "none"
plt.rcParams["svg.fonttype"] = "none"

os.chdir("/data/yosef2/Th1_IL23/analysis/cd160_knockout")


# %% Load the data

tpm = pd.read_table("rnaseq/tpm.txt", index_col=0)
meta = pd.read_table("rnaseq/meta.txt", index_col=0)

# %% Drop S1 - outlier sample

valid_samples = meta.drop("S1", axis=0).index
tpm = tpm.loc[:, valid_samples]
meta = meta.loc[valid_samples]


# %% Load DE Results

de_results = pd.read_table("de_genotype.txt", index_col=0)
de_results = de_results.drop("TDTOMATO", axis=0)

# %% Load cluster 1vall results
cluster_1vall_results_file = "../10xAll/combined_LPL/cluster_together/cluster_de_1vAll_edgeR//cluster_de_results.txt.gz"
cluster_1vall_results = pd.read_table(cluster_1vall_results_file)


# %% Compute per-cluster signatures


def get_sigs(df):
    sigs = {}
    for cluster in df["cluster"].unique():
        sub = df.loc[
            lambda x: (
                (x["cluster"] == cluster)
                & (x["FDR"] < 0.05)
                & (x["logFC"].abs() > 0.25)
            )
        ]

        sig_up_genes = sub.loc[lambda x: x["logFC"] > 0].GeneSymbol.tolist()
        sig_dn_genes = sub.loc[lambda x: x["logFC"] < 0].GeneSymbol.tolist()

        sig_up_genes = [x.upper() for x in sig_up_genes]
        sig_dn_genes = [x.upper() for x in sig_dn_genes]

        sig = {x: 1 for x in sig_up_genes}
        sig.update({x: -1 for x in sig_dn_genes})
        sig = pd.Series(sig)
        sigs[cluster] = sig

    return sigs


sigs_1vall = get_sigs(cluster_1vall_results)
sigs_1vall = {k + "-1vall": v for k, v in sigs_1vall.items()}


# Score these signatures on the cd160 samples

log_tpm = np.log(tpm + 1)
log_tpm_norm = log_tpm.subtract(log_tpm.mean(axis=0), axis=1).divide(
    log_tpm.std(axis=0), axis=1
)
log_tpm_norm = log_tpm_norm.drop(
    "CD160", axis=0
)  # Can't use this to test as we knocked it out specifically


def score_sig(sig_name, sig_values):
    sig_up_genes = sig_values.index[sig_values > 0] & log_tpm_norm.index
    sig_dn_genes = sig_values.index[sig_values < 0] & log_tpm_norm.index

    sig_up_score = log_tpm_norm.loc[sig_up_genes].mean(axis=0)
    sig_dn_score = log_tpm_norm.loc[sig_dn_genes].mean(axis=0)
    sig_score = sig_up_score - sig_dn_score
    sig_score.name = sig_name
    return sig_score


sig_scores = [score_sig(name, vals) for name, vals in sigs_1vall.items()]
sig_scores = pd.concat(sig_scores, axis=1)

plot_data = sig_scores.join(meta)

from scipy.stats import ttest_ind


pvals = []
sigs_to_test = ["c2-1vall", "c7-1vall", "c9-1vall"]
for sig in sigs_to_test:
    a = plot_data.loc[lambda x: x["Genotype"] == "ctrl", sig]
    b = plot_data.loc[lambda x: x["Genotype"] == "KO", sig]

    ttest_result = ttest_ind(a.values, b.values)
    pvals.append(ttest_result.pvalue)

from statsmodels.stats.multitest import multipletests

pvals_corrected = multipletests(pvals, method="fdr_bh")[1]
pvals_corrected = pd.Series(pvals_corrected, index=sigs_to_test)


fig, axs = plt.subplots(1, 3, figsize=(7, 5))

genotype_map = {
    'ctrl': 'WT',
    'KO': '$CD160^{-/-}$',
}
plot_data['GenotypeNice'] = plot_data['Genotype'].map(genotype_map)
for sig, ax in zip(["C2", "C7", "C9"], axs.ravel()):
    plt.sca(ax)
    sns.boxplot(
        data=plot_data,
        x="GenotypeNice",
        y=sig.lower() + "-1vall",
        palette=["#AAAAAA", colors[0]],
        whis=99,
        width=0.5,
    )
    plt.xlabel("")
    plt.ylabel("")
    plt.title(sig + " vs. All")
    if ax == axs[0]:
        plt.ylabel("Signature Score")

plt.subplots_adjust(wspace=0.55)
plt.savefig("Figures/sig_scores_1vall.svg")


# # Now do this in reverse


# Define the cd160 signature

up_genes = de_results.loc[
    lambda x: ((x["FDR"] < 0.1) & (x["logFC"] > 0.5))
].index

dn_genes = de_results.loc[
    lambda x: ((x["FDR"] < 0.1) & (x["logFC"] < 0.5))
].index

up_genes = up_genes.str.capitalize()
dn_genes = dn_genes.str.capitalize()

print(len(up_genes), len(dn_genes))


# Load single-cell gexp and define signature scores

scaled_exp = feather.read_dataframe(
    "../10xAll/combined_LPL/data/expression_scaled.feather"
)
scaled_exp.index = scaled_exp["index"]
scaled_exp.index = [x.capitalize() for x in scaled_exp.index]
scaled_exp = scaled_exp.drop("index", axis=1)

clusters = pd.read_table(
    "../10xAll/combined_LPL/cluster_together/clusters.txt", index_col=0
)
clusters = clusters.loc[scaled_exp.columns]
cluster_colors = json.load(
    open("../10xAll/combined_LPL/cluster_together/cluster_colors.json")
)

umap = pd.read_table("../10xAll/combined_LPL/umap/umap.txt", index_col=0)
umap = umap.loc[scaled_exp.columns]

# %% Evaluate a signature score

logexp = np.log(scaled_exp + 1)

up_score = logexp.loc[[x in up_genes for x in logexp.index], :].mean(axis=0)
dn_score = logexp.loc[[x in dn_genes for x in logexp.index], :].mean(axis=0)

sig_score_lpl = up_score - dn_score


scores_df = clusters.join(sig_score_lpl.rename("SigScore"))

scores_df["cluster_group_c9"] = "Other"
scores_df.loc[lambda x: x["Cluster"] == "c9", "cluster_group_c9"] = "c9"

scores_df["cluster_group_c2"] = "Other"
scores_df.loc[lambda x: x["Cluster"] == "c2", "cluster_group_c2"] = "c2"

scores_df["cluster_group_c7"] = "Other"
scores_df.loc[lambda x: x["Cluster"] == "c7", "cluster_group_c7"] = "c7"


cluster_order = (
    scores_df.groupby("Cluster")["SigScore"].median().sort_values().index
)

plt.figure(figsize=(7, 4))

sns.boxplot(
    x="Cluster",
    y="SigScore",
    data=scores_df,
    order=cluster_order,
    width=0.65,
    whis=99,
)
plt.ylabel("$CD160^{-/-}$ Signature Score")
plt.savefig("Figures/sig_scores_cd160ko_percluster.svg")


fig, axs = plt.subplots(1, 3, figsize=(7, 4))

for sig, ax in zip(["c2", "c7", "c9"], axs.ravel()):
    plt.sca(ax)
    sns.boxplot(
        data=scores_df,
        x="cluster_group_" + sig,
        y="SigScore",
        palette=["#AAAAAA", colors[0]],
        whis=99,
        width=0.5,
    )
    plt.xlabel("")
    plt.ylabel("")

    scores_num = scores_df.loc[clusters.Cluster == sig]["SigScore"]
    scores_other = scores_df.loc[clusters.Cluster != sig]["SigScore"]

    test_result = ttest_ind(scores_num, scores_other)
    print(sig, test_result)

    if ax == axs[0]:
        plt.ylabel("$CD160^{-/-}$ Signature Score")

plt.subplots_adjust(wspace=0.55)
plt.savefig("Figures/sig_scores_cd160ko_1vall.svg")
