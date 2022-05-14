import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams['svg.fonttype'] = 'none'

# %% Load the differential expression sets


de_th17_ctrl = (
    pd.read_table("../Th17_Pipeline/DE_Th17_ctrl_noqc/results.txt.gz")
    .loc[lambda x: x['contrast'] == 'GFPeGFP_pos']
    .set_index('primerid')
)

de_th17_ko = (
    pd.read_table("../Th17_Pipeline/DE_Th17_KO_noqc/results.txt.gz")
    .loc[lambda x: x['contrast'] == 'GFPeGFP_pos']
    .set_index('primerid')
)

de_th1_ctrl = (
    pd.read_table("../Th1_Pipeline/DE_Th1_ctrl_noqc/results.txt.gz")
    .loc[lambda x: x['contrast'] == 'GFPeGFP_pos']
    .set_index('primerid')
)

de_th1_ko = (
    pd.read_table("../Th1_Pipeline/DE_Th1_KO_noqc/results.txt.gz")
    .loc[lambda x: x['contrast'] == 'GFPeGFP_pos']
    .set_index('primerid')
)

# Fix the NaNs in the logFC (for sorting and grouping) - just get the sign correct

for x in [de_th17_ctrl, de_th17_ko, de_th1_ctrl, de_th1_ko]:
    to_fix = pd.isnull(x['logFC'])
    x.loc[to_fix, 'logFC'] = x.loc[to_fix, 'coefD']

# %% Load in the gene expression

tpm = pd.read_table(
    "../data/tpm.txt.gz", index_col=0
)

meta17 = pd.read_table(
    "../Th17_Pipeline/data_filtered/meta.txt.gz", index_col=0
)

meta1 = pd.read_table(
    "../Th1_Pipeline/data_filtered/meta.txt.gz", index_col=0
)

meta = pd.concat((meta1, meta17), axis=0)
tpm = tpm.loc[:, meta.index]

# %% Organize the gene expression into the appropriate groupings

# Maps (Stimulus, Genotype, GFP) to a data matrix - already z-normalized within the Stim/Geno grouping
tpm_norm_all = {}  

for stim in meta.Stimulus.unique():
    for geno in meta.Genotype.unique():

        data = tpm.loc[
            :, (
                (meta['Stimulus'] == stim) &
                (meta['Genotype'] == geno)
            )]

        data = data.subtract(
            data.mean(axis=1), axis=0
        ).divide(
            data.std(axis=1), axis=0
        )

        data = data.fillna(0)

        for gfp in meta.GFP.unique():

            data_sub = data.loc[
                :, (meta.loc[data.columns, 'GFP'] == gfp)
            ]

            data_sub = data_sub.sample(frac=1, axis=1)  # Randomize the ordering

            tpm_norm_all[(stim, geno, gfp)] = data_sub


# %% Divide genes into groups based on test results


# Common effect

both = (
    (
        de_th17_ctrl.index[de_th17_ctrl.FDR < 0.05] &
        de_th1_ctrl.index[de_th1_ctrl.pvalH < 0.05]
    ) | (
        de_th1_ctrl.index[de_th1_ctrl.FDR < 0.05] &
        de_th17_ctrl.index[de_th17_ctrl.pvalH < 0.05]
    )
)

both_common = both[
    (
        de_th1_ctrl.loc[both].logFC *
        de_th17_ctrl.loc[both].logFC
    ) > 0
]

both_divergent = both.difference(both_common)


th1_only = de_th1_ctrl.index[de_th1_ctrl.FDR < 0.05].difference(both)
th17_only = de_th17_ctrl.index[de_th17_ctrl.FDR < 0.05].difference(both)


# %% Plot a group of genes
def plot_group(genes_to_plot, ticklabel_size=5, genes_to_label=None, label_position='left', inches_per_gene=.2, fig_width=7, fig_v_margin=1):

    cmap = LinearSegmentedColormap.from_list(
        'purple_black_yellow', ['#FF00FF', '#000000', '#FFFF00']
    )

    groups = [
        ('12+21+23', 'ctrl'),
        ('1+6+23', 'ctrl'),
    ]

    if genes_to_label is None:
        genes_to_label == set(genes_to_plot)

    fig_height = len(genes_to_plot) * inches_per_gene + fig_v_margin

    fig = plt.figure(figsize=(fig_width, fig_height))

    gs0 = fig.add_gridspec(1, 2, wspace=.3)

    i = 0
    for gs, (stim, geno) in zip(gs0, groups):

        gs_sub = gs.subgridspec(1, 2, wspace=.1)
        ax_neg = fig.add_subplot(gs_sub[0])
        ax_pos = fig.add_subplot(gs_sub[1])

        data_neg = tpm_norm_all[(stim, geno, 'eGFP_neg')]
        data_pos = tpm_norm_all[(stim, geno, 'eGFP_pos')]

        sns.heatmap(
            data_neg.loc[genes_to_plot], xticklabels=False, yticklabels=False,
            ax=ax_neg, cmap=cmap,
            vmin=-1.5, vmax=1.5, cbar=False, rasterized=True
        )
        ax_neg.set_xlabel('')
        ax_neg.set_ylabel('')
        ax_neg.set_title('{}_{}\n{}'.format(stim, geno, 'GFP-'), size=8)

        if i == 0 and label_position == 'left':
            ticks_and_labels = [
                (i, gene) for i, gene in enumerate(genes_to_plot) if gene in genes_to_label
            ]
            ticks = [x[0] for x in ticks_and_labels]
            labels = [x[1] for x in ticks_and_labels]
            ax_neg.set_yticks(ticks)
            ax_neg.set_yticklabels(labels, size=ticklabel_size)

        sns.heatmap(
            data_pos.loc[genes_to_plot], xticklabels=False, yticklabels=False,
            ax=ax_pos, cmap=cmap,
            vmin=-1.5, vmax=1.5, cbar=False, rasterized=True
        )
        ax_pos.set_xlabel('')
        ax_pos.set_ylabel('')
        ax_pos.set_title('{}_{}\n{}'.format(stim, geno, 'GFP+'), size=8)

        if i == 1 and label_position == 'right':
            ticks_and_labels = [
                (i, gene) for i, gene in enumerate(genes_to_plot) if gene in genes_to_label
            ]
            ticks = [x[0] for x in ticks_and_labels]
            labels = [x[1] for x in ticks_and_labels]
            ax_pos.yaxis.set_ticks_position("right")
            ax_pos.set_yticks(ticks)
            ax_pos.set_yticklabels(labels, size=ticklabel_size)

        i = i + 1

    plt.subplots_adjust(bottom=fig_v_margin / fig_height, top=1 - fig_v_margin / fig_height)

    return fig


# %%

# Change this to None to label all genes
genes_to_label = {
    "IL23R", "EGFP_IRES_IL23R", "IL22", "ERMN", "PERP", "S1PR1",
    "GPR18", "PDE4B", "SMOX", "LMNA", "STOM", "KLF6", "RORC",
    "SOCS3", "RORA", "DUSP3", "MT2", "FLOT1", "ANXA2", "BATF",
    "CASP6", "STAT5A", "PTPN1", "S100A11", "CD81", "IL12RB2", "IL1R1",
    "TNFRSF18", "TNFRSF4", "LTB", "IFI27L2A", "SLAMF6", "LAMP1",
    "TRIM27", "CD3E", "MAP3K1", "IL27RA", "IL16", "GATA3", "IL4",
    "IL17F", "IL17A", "BHLHE40", "JUNB", "CCR3", "DUSP2", "IFNG",
    "IL18R1", "TUSC2", "PELI1", "PLEK", "MAPK1", "CSF1", "FURIN",
    "MYC", "MAX", "IL2RA", "CD72",
}

inches_per_gene = .02  # vertical scaling
ticklabel_size = 8
fig_width = 4
fig_v_margin = .5

gene_lists = {}  # Saving these here so that the ordering is preserved as well

genes_to_plot = (
    de_th1_ctrl.loc[th1_only].sort_values("logFC", ascending=False).index
)
fig = plot_group(
    genes_to_plot, ticklabel_size, genes_to_label=genes_to_label,
    label_position="right", inches_per_gene=inches_per_gene,
    fig_width=fig_width, fig_v_margin=fig_v_margin
)
fig.savefig("heatmap_th1_only.svg", dpi=300)
gene_lists["th1_only"] = genes_to_plot

genes_to_plot = (
    de_th17_ctrl.loc[th17_only].sort_values("logFC", ascending=False).index
)
fig = plot_group(
    genes_to_plot, ticklabel_size, genes_to_label=genes_to_label,
    label_position="right", inches_per_gene=inches_per_gene,
    fig_width=fig_width, fig_v_margin=fig_v_margin
)
fig.savefig("heatmap_th17_only.svg", dpi=300)
gene_lists["th17_only"] = genes_to_plot

genes_to_plot = (
    de_th1_ctrl.loc[both_common].sort_values("logFC", ascending=False).index
)
fig = plot_group(
    genes_to_plot, ticklabel_size, genes_to_label=genes_to_label,
    label_position="left", inches_per_gene=inches_per_gene,
    fig_width=fig_width, fig_v_margin=fig_v_margin
)
fig.savefig("heatmap_both_common.svg", dpi=300)
gene_lists["both_common"] = genes_to_plot

genes_to_plot = (
    de_th1_ctrl.loc[both_divergent].sort_values("logFC", ascending=False).index
)
fig = plot_group(
    genes_to_plot, ticklabel_size, genes_to_label=genes_to_label,
    label_position="left", inches_per_gene=inches_per_gene,
    fig_width=fig_width, fig_v_margin=fig_v_margin
)
fig.savefig("heatmap_both_divergent.svg", dpi=300)
gene_lists["both_divergent"] = genes_to_plot

# Save sheets

writer = pd.ExcelWriter('heatmap_gene_lists.xlsx')

for listname, genes in gene_lists.items():
    pd.DataFrame({'Genes': genes}) \
        .to_excel(writer, sheet_name=listname, header=False, index=False)

writer.save()
