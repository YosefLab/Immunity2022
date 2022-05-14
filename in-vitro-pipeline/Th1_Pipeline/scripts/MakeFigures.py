import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib import gridspec, cm
import numpy as np
import pandas as pd
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

# %% Make the Heatmaps

tpm = pd.read_table("normalized_noqpc/tpm.txt.gz", index_col=0)
meta = pd.read_table("data_filtered/meta.txt.gz", index_col=0)

tpm = tpm.loc[
    :, meta.index[meta.Genotype == 'ctrl']
]

meta = meta.loc[meta.Genotype == 'ctrl']

plot_data = np.log2(tpm+1)
plot_data = plot_data.subtract(plot_data.mean(axis=1), axis=0) \
                     .divide(plot_data.std(axis=1), axis=0)

plot_data_pos = plot_data.loc[
    :, meta.index[meta.GFP == 'eGFP_pos']
]

plot_data_neg = plot_data.loc[
    :, meta.index[meta.GFP == 'eGFP_neg']
]

# %%


genes = [
    ['IL23R', 'IFNG', 'CSF1', 'IL22'],
    ['IL18R1', 'IL2RA', 'IL1R2', 'SOCS3'],
    ['CD3E', 'LGALS1', 'CD81', 'S1PR1'],
    ['STAT5A', 'BATF', 'RORC', 'RORA'],
    ['IL4', 'CCR8', 'GATA3']
]

cmap = LinearSegmentedColormap.from_list("purple", ["#FF00FF", "#000000", "#FFFF00"])

fig = plt.figure(figsize=(5, 7))
gs0 = gridspec.GridSpec(
    len(genes), 2, height_ratios=[len(x) for x in genes], wspace=.05
)

for i, gsub in enumerate(genes):

    ax = plt.subplot(gs0[i, 0])

    sns.heatmap(
        plot_data_neg.loc[gsub], cmap=cmap,
        vmin=-2, vmax=2, yticklabels=True, xticklabels=False,
        ax=ax, cbar=False, rasterized=True
    )

    ax = plt.subplot(gs0[i, 1])

    sns.heatmap(
        plot_data_pos.loc[gsub], cmap=cmap,
        vmin=-2, vmax=2, yticklabels=False, xticklabels=False,
        ax=ax, cbar=False, rasterized=True
    )

plt.subplots_adjust(left=.15, right=0.8)

cax = fig.add_axes([.55, .05, .2, .02])
norm = Normalize(vmin=-2, vmax=2)
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
             orientation='horizontal')

fig.patch.set_visible(False)
plt.savefig('Figures/GFP_DE_Heatmap.svg', dpi=300)

# %% Make the Heatmaps - Ctrl AND KO

tpm = pd.read_table("normalized_noqpc/tpm.txt.gz", index_col=0)
meta = pd.read_table("data_filtered/meta.txt.gz", index_col=0)

# Regular normalization
# plot_data = np.log2(tpm+1)
# plot_data = plot_data.subtract(plot_data.mean(axis=1), axis=0) \
#                      .divide(plot_data.std(axis=1), axis=0)

# Genotype-specific normalization
tpm_c = tpm.loc[
    :, meta.index[meta.Genotype == 'ctrl']
]

plot_data_c = np.log2(tpm_c+1)
plot_data_c = plot_data_c.subtract(plot_data_c.mean(axis=1), axis=0) \
                     .divide(plot_data_c.std(axis=1), axis=0)

tpm_k = tpm.loc[
    :, meta.index[meta.Genotype == 'KO']
]

plot_data_k = np.log2(tpm_k+1)
plot_data_k = plot_data_k.subtract(plot_data_k.mean(axis=1), axis=0) \
                     .divide(plot_data_k.std(axis=1), axis=0)

plot_data = pd.concat((plot_data_c, plot_data_k), axis=1)


# %%

genes = [
    ['IL22', 'GZMB', 'CCR3', 'LAMP1',
     'ITGB1', 'GPR18', 'RAC1', 'PELI1',
     'MAPK1', 'MAP2K2', 'DUSP2', 'TRIM27',
     'JUND', 'LTA', 'BHLHE40', 'RBPJ']
]

columns = [
    ('ctrl', 'eGFP_neg'),
    ('ctrl', 'eGFP_pos'),
    ('KO', 'eGFP_neg'),
    ('KO', 'eGFP_pos'),
]

cmap = LinearSegmentedColormap.from_list("purple", ["#FF00FF", "#000000", "#FFFF00"])

fig = plt.figure(figsize=(10, 7))
gs0 = gridspec.GridSpec(
    len(genes), 4, height_ratios=[len(x) for x in genes], wspace=.05
)

for i, gsub in enumerate(genes):

    for j, (geno, gfp) in enumerate(columns):

        ax = plt.subplot(gs0[i, j])

        plot_data_sub = plot_data.loc[
            :, meta.index[
                (meta.GFP == gfp) &
                (meta.Genotype == geno)
            ]
        ]

        sns.heatmap(
            plot_data_sub.loc[gsub], cmap=cmap,
            vmin=-2, vmax=2, yticklabels=True, xticklabels=False,
            ax=ax, cbar=False, rasterized=True
        )

        if j > 0:
            plt.yticks([])

        plt.title('{} - {}'.format(geno, gfp))


plt.subplots_adjust(left=.15, right=0.8)

cax = fig.add_axes([.55, .05, .2, .02])
norm = Normalize(vmin=-2, vmax=2)
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
             orientation='horizontal')

fig.patch.set_visible(False)
plt.show()
#plt.savefig('Figures/GFP_DE_Heatmap.svg', dpi=300)
