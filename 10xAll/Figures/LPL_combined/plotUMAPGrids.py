import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import feather
import common
from itertools import zip_longest

plt.rcParams['svg.fonttype'] = 'none'

# %% Load Data

expression = feather.read_dataframe(
    "../../combined_LPL/data/expression_scaled.feather"
).set_index("index")/3000*10000

meta = pd.read_table("../../combined_LPL/data/meta.txt.gz", index_col=0)
umap = pd.read_table("../../combined_LPL/umap/umap.txt", index_col=0)
clusters = pd.read_table(
    "../../combined_LPL/cluster_together/clusters.txt", index_col=0
)

meta = meta.join(umap).join(clusters)

# %% Plot 3x2 grid of cell-cycle genes

fig, axs = plt.subplots(
    2, 3, gridspec_kw=dict(left=0.1, right=0.85, top=0.9, bottom=0.1),
    figsize=(8.5*1.2, 5*1.2)
)


cbar_ax = fig.add_axes([0.68, 0.1, 0.01, 0.1])

def despine(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.sca(ax)
    plt.xticks([])
    plt.yticks([])


genes = ["LPL-8", "Mki67", "Cdc6", "Cdk1", "Cdc20", "Ccnb1"]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="grays", colors=["#dddddd", "000000"]
)

kwargs = dict(cmap=cmap, vmin=0, vmax=5, s=1.5, edgecolors="none")

for g, ax in zip_longest(genes, axs.ravel()):

    if g is None:
        ax.remove()
        continue

    plt.sca(ax)

    if g == "LPL-8":

        c_color = common.colors.cmap_lpl_clusters["c8"]

        c = [c_color if x == "c8" else "#dddddd" for x in meta.Cluster]

        plt.scatter(x=meta.umap1, y=meta.umap2, c=c, s=1.5, rasterized=True, edgecolors="none")
        plt.title(common.name_map_lpl_clusters["c8"], size=13)

    else:

        exp = np.log2(expression.loc[g] + 1)
        exp = exp.loc[meta.index]

        pc = plt.scatter(
            x=meta.umap1, y=meta.umap2, c=exp, **kwargs, rasterized=True
        )

    despine(ax)
    plt.title(g, size=11)

plt.colorbar(
    pc, cax=cbar_ax, ticks=[0, 2.5, 5], label="Expression\n$log_2(CP10K+1)$"
)

# plt.show()

plt.savefig('umap_grid_cellcycle.svg', dpi=600)


# %% Plot 4x2 grid of expected/control surface markers

fig, axs = plt.subplots(
    2, 4, gridspec_kw=dict(left=0.08, right=0.88, top=0.9, bottom=0.1),
    figsize=(10.625*1.2, 5*1.2)
)


cbar_ax = fig.add_axes([0.90, 0.1, 0.01, 0.1])

def despine(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.sca(ax)
    plt.xticks([])
    plt.yticks([])


genes = [
    "Cd3e",
    "Ptprc",
    "Ncam1",
    "Fcgr3",
    "Cd3g",
    "Cd4",
    "Klrb1",
    "Cd8a",
]

cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="grays", colors=["#dddddd", "000000"]
)

kwargs = dict(cmap=cmap, vmin=0, vmax=5, s=1.5, edgecolors="none")

for g, ax in zip_longest(genes, axs.ravel()):

    if g is None:
        ax.remove()
        continue

    plt.sca(ax)

    exp = np.log2(expression.loc[g] + 1)
    exp = exp.loc[meta.index]

    pc = plt.scatter(
        x=meta.umap1, y=meta.umap2, c=exp, **kwargs, rasterized=True
    )

    despine(ax)
    plt.title(g, size=11)

plt.colorbar(
    pc, cax=cbar_ax, ticks=[0, 2.5, 5], label="Expression\n$log_2(CP10K+1)$"
)

# plt.show()

plt.savefig('umap_grid_celltype_markers.svg', dpi=600)
