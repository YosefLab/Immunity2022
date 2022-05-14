import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

plt.rcParams['svg.fonttype'] = 'none'

# %% Read data

df = pd.read_table("../../Targets/Targets.txt.gz", index_col=0)

# %% Compute Summary columns for figure

# In-vivo cluster specific

a1 = df[[
    "c2_Sig",
    "c2_Geno_Sig"
]]

a2 = df[[
    "c9_Sig",
    "c9_Geno_Sig"
]]/2

in_vivo_1 = pd.concat((a1, a2), axis=1).max(axis=1).clip(upper=2)

# In-vivo tissue markers

def resolve_tissue_de(logfc, fdr):
    if fdr < 0.1 and logfc > 0:
        return logfc
    return 0.0


# iel_tissue_marker = {
#     gene: resolve_tissue_de(logfc, fdr)
#     for gene, logfc, fdr in zip(df.index, df["logFC_IEL_TissueDE"], df["FDR_IEL_TissueDE"])
# }
#
# iel_tissue_marker = pd.Series(iel_tissue_marker)

lpl_tissue_marker = {
    gene: resolve_tissue_de(logfc, fdr)
    for gene, logfc, fdr in zip(df.index, df["logFC_LPL_TissueDE"], df["FDR_LPL_TissueDE"])
}

lpl_tissue_marker = pd.Series(lpl_tissue_marker)

a3 = df[[
    "LPL_Geno_Sig",
    # "IEL_Geno_Sig",
]]

in_vivo_2 = pd.concat((lpl_tissue_marker,  a3), axis=1) \
    .max(axis=1).clip(upper=2)


df["In-Vivo Cluster"] = in_vivo_1
df["In-Vivo Tissue"] = in_vivo_2

# In-Vitro

df["In-Vitro Summary"] = (
    df[["In-vitro Th1 Dep Sig", "In-vitro Th1"]].max(axis=1).clip(lower=0, upper=3)
    #+ (df["IL23R_KO_Th17_Microarray"]*-1).clip(lower=0, upper=2)/5
)

# Zero out IL23R in the in-vitro tests as the GFP sorting is designed
# to capture this transcript in differential comparisons
df.loc["Il23r", "In-Vitro Summary"] = 0

# df["In-Vitro Summary"] = (
#     df[["In-vitro Th1 Dep Sig"]].max(axis=1).clip(lower=0, upper=2) +
#     df["In-vitro_Th1_IL23R-dep_TF"].clip(upper=2, lower=0)/5 +
#     (df["IL23R_KO_Th17_Microarray"]*-1).clip(lower=0, upper=2)/5
# )/2


df["Score"] = (
    df["In-Vivo Cluster"]*2 + df["In-Vivo Tissue"]*1 +
    df["In-Vitro Summary"]*1 + df["GWAS?"]*2
)

df = df.sort_values("Score", ascending=False)

# %% Plot the things


def plot_v3():

    N = 40

    fig = plt.figure(figsize=(15, 3.5))

    gs0 = fig.add_gridspec(2, 1, height_ratios=[1.6, 1], hspace=.4)
    gsTop = gs0[0].subgridspec(4, 1, hspace=.5)

    cbar_ax = fig.add_axes(
        [.6, .905, .1, .02]
    )

    rows = {
        "In-Vivo Cluster-Specific": df["In-Vivo Cluster"],
        "In-Vivo Tissue-Specific": df["In-Vivo Tissue"],
        "In-Vitro Summary": df["In-Vitro Summary"],
        "GWAS": df["GWAS?"],
    }

    cmap_rest = "Greens"
    cmap_gwas = LinearSegmentedColormap.from_list(
        "gwas", ["#FFFFFF", "#9e62a6"]
    )

    for i, row in enumerate(rows.keys()):

        ax = fig.add_subplot(gsTop[i, 0])
        plt.sca(ax)

        plot_data = rows[row].iloc[0:N].rename(row).to_frame().T
        sns.heatmap(
            plot_data,
            cbar=False,
            xticklabels=True,
            yticklabels=True,
            cmap=cmap_rest if row != "GWAS" else cmap_gwas,
            vmin=0,
            vmax=plot_data.values.max(),
            linewidths=.5,
            linecolor='#CCCCCC',
        )

        ax.tick_params(
            axis="y", left=False, right=True, labelleft=False, labelright=True
        )
        ax.tick_params(
            axis="x", top=False, bottom=(i == len(rows)-1), labeltop=False, labelbottom=False
        )
        plt.yticks(rotation=0)

    plt.colorbar(
        ScalarMappable(norm=Normalize(0, 1), cmap=cmap_rest),
        cax=cbar_ax,
        orientation='horizontal',
        ticks=[0, 1],
    )
    cbar_ax.set_xticklabels(['Min', 'Max'])
    cbar_ax.tick_params(
        axis='x',
        bottom=False, labelbottom=False,
        top=True, labeltop=True,
    )

    ax = fig.add_subplot(gs0[1, 0])
    plt.sca(ax)

    plt.bar(np.arange(N) + 0.5, df["Score"][0:N].values,
            tick_label=df.index[0:N])
    plt.ylabel("Score")
    plt.ylim(0, df["Score"].max() * 1.05)
    plt.xticks(rotation=90, size=8)

    # plt.yticks([1.5, 2, 2.5])

    plt.subplots_adjust(right=0.7, bottom=0.3, hspace=0.6, top=.9)

    for x in ax.get_xticklabels():
        if x.get_text() in {'Cd160', 'S1pr1', 'Cd226', 'Zfp36l1', 'Gpr18'}:
            x.set_color('#0b4967')

    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_facecolor('#F4F4F4')
    ax.set_axisbelow(True)
    plt.grid(axis='y', color='#FFFFFF')
    plt.xlim(0, N)

    ax.tick_params(
        axis="y", left=False, right=True, labelleft=False, labelright=True
    )
    ax.yaxis.set_label_position("right")

    return fig

plot_v3()
plt.savefig('targets.svg')


# %% Also, create the supplementary table that goes with this figure

df_out = df[
    [
        "c9_logFC",
        "c9_FDR",
        "c9_Geno_logFC",
        "c9_Geno_FDR",
        "c2_logFC",
        "c2_FDR",
        "c2_Geno_logFC",
        "c2_Geno_FDR",
        "logFC_LPL_TissueDE",
        "FDR_LPL_TissueDE",
        "LPL_Geno_logFC",
        "LPL_Geno_FDR",
        "In-Vivo Cluster",
        "In-Vivo Tissue",
        "In-Vitro Summary",
        "GWAS?",
        "Score",
    ]
].copy()

df_out.to_csv("targets_supplementary_table_pre.txt", sep="\t")
