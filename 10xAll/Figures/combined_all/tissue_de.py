import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import feather
import common
import bio_utils.plots

plt.rcParams['svg.fonttype'] = 'none'

# %% Load counts to get averages...

scaled = feather.read_dataframe("analysis/data/expression_scaled.feather").set_index('index')
meta = pd.read_table("analysis/data/meta.txt.gz", index_col=0)
gene_info = pd.read_table("analysis/data/genes.tsv", index_col=0)

tissue_means = scaled.T.join(meta.Tissue).groupby('Tissue').mean().T / 3000 * 10000


# %% Load the cpm
cpm = pd.read_table("analysis/TissueDE/tissue_de2_cpm.txt.gz", index_col=0)
cpm = cpm.groupby('GeneSymbol').sum()
cpm = cpm.loc[tissue_means.index]
cpm = cpm.loc[:, tissue_means.columns]

# %% Load top DE genes
de_all_res = pd.read_table("analysis/TissueDE/tissue_de2.txt.gz", index_col=0)

genes = pd.Index(de_all_res.sort_values('LR', ascending=False).GeneSymbol[0:500])

data = cpm.loc[genes]
data = data.divide(data.max(axis=1), axis=0)

# %% plot specific genes

genes = pd.Index(de_all_res.sort_values('LR', ascending=False).GeneSymbol[0:500])

data = cpm.loc[genes]
data = data.divide(data.max(axis=1), axis=0)

genes_to_show = [
    #"Vps37b",
    #"Ubald2",
    "Rgs1",
    #"Tnfaip3",
    #"Junb",
    #"Bhlhe40",
    "Ctla4",
    #"Jund",
    #"Ifngr1",
    "Tnfrsf4",
    #"Crem",
    #"Dusp5",
    #"Icos",
    "Itgb1",
    #"Dusp1",
    "Cd52",
    #"Ramp3",
    "Cd48",
    #"Pde4b",
    #"Nfkbiz",
    #"Kdm6b",
    #"Fth1",
    #"Fosl2",
    "Lgals1",
    "Klf2",
    #"Rasgrp2",
    "Itgb2",
    "Cd69",
    "Itgb7",
    "S100a4",
    #"Dnajb1",
    #"Rora",
    #"Tnfrsf18",
    #"Fos",
    "Il4ra",
    "Cxcr6",
    #"Klf6",
    "Cd28",
    #"Zfp36l2",
    "Ly6c2",
    "Fasl",
    "Itgae",
    "Gpr183",
    #"Fosb",
    #"Tnfsf8",
    #"Jun",
    "Ifnar1",
    "Tnfsf11",
    "Cd160",
    #"Anxa1",
    #"Plac8",
    #"S1pr4",
    #"Ifng",
    #"Dusp2",
    "Cxcl10",
    #"Zfp36l1",
    #"Stat3",
    "Ccr9",
    "Il2rg",
    #"Il22",
    #"Gzmb",
    "Ccr5",
    #"Tnfrsf9",
    #"Dusp4",
    #"Socs3",
    "Cd7",
    #"Ahr",
    "Lag3",
    #"Klf3",
    "Cxcr3",
    #"Gzma",
    "S1pr1",
    "Cd53",
    "Il12rb2",
    "Gpr18",
]

from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.colors import LinearSegmentedColormap
Z = linkage(data, metric='cosine', method='average')
ii = leaves_list(Z)

# %% Plot heatmap

fig = plt.figure(figsize=(4, 8))
cbax = fig.add_axes([.05, .1, .02, .2])
ax = fig.add_axes([.2, .1, .5, .8])
plt.sca(ax)
data = data.loc[:, common.order_tissue]

cmap = LinearSegmentedColormap.from_list(
    name="PuYl", colors=["#DD00DD", "#000000", "#DDDD00"])

sns.heatmap(
    data.iloc[ii],
    vmin=0, vmax=1, cmap=cmap,
    yticklabels=True,
    xticklabels=[common.name_map_tissue[x] for x in data.columns],
    cbar_ax=cbax, rasterized=True,
)

ticks = []
labels = []
ax = plt.gca()
for tick, tlb in zip(ax.get_yticks(), ax.get_yticklabels()):
    if tlb.get_text() in genes_to_show:
        labels.append(tlb.get_text())
        ticks.append(tick)

plt.yticks(ticks, labels, rotation=0)
plt.xlabel("")
plt.ylabel("")
plt.xticks(rotation=45)
plt.tick_params(axis='y', 
                labelleft=False, labelright=True,
                left=False, right=True,
                )

plt.savefig('de_tissue_pre.svg', dpi=300)

# %% Or, try dot again but with scaled

d_long = data.loc[genes_to_show] \
    .reset_index() \
    .melt(id_vars='GeneSymbol', value_name='Expression')
y_order = [x for x in data.index[ii] if x in genes_to_show]

plt.figure(figsize=(5, 9))
bio_utils.plots.dot_plot(
    data=d_long, y='GeneSymbol', x='Tissue', size='Expression',
    y_order=y_order, x_order=common.order_tissue,
    lpad=0.2, size_norm=(0, .9),
)
plt.xticks(rotation=90)
# plt.show()
plt.savefig('de_tissue_dot_pre.svg')
