import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import feather
import common

plt.rcParams['svg.fonttype'] = 'none'

# %% Load Data

expression = feather.read_dataframe(
    "analysis_lpl/data/expression_scaled.feather"
).set_index("index")/3000*10000

meta = pd.read_table("analysis_lpl/data/meta.txt.gz", index_col=0)
umap = pd.read_table("analysis_lpl/umap/umap.txt", index_col=0)
clusters = pd.read_table(
    "analysis_lpl/cluster_together/clusters.txt", index_col=0
)

meta = meta.join(umap).join(clusters)

# %% Plot the grids!

fig = plt.figure(figsize=(7, 6))

gs = fig.add_gridspec(3, 3, left=0.1, right=0.85, top=0.9, bottom=0.1)

ax1 = fig.add_subplot(gs[0:2, 0:2])

ax1s = fig.add_subplot(gs[0, 2])
ax2s = fig.add_subplot(gs[1, 2])
ax3s = fig.add_subplot(gs[2, 2])
ax4s = fig.add_subplot(gs[2, 1])
ax5s = fig.add_subplot(gs[2, 0])

cbar_ax = fig.add_axes([0.88, 0.1, 0.01, 0.1])

small_axs = [ax1s, ax2s, ax3s, ax4s, ax5s]


def despine(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.sca(ax)
    plt.xticks([])
    plt.yticks([])


# Plot the big umap

plt.sca(ax1)

c_color = common.colors.cmap_lpl_clusters["c7"]

c = [c_color if x == "c7" else "#dddddd" for x in meta.Cluster]

plt.scatter(x=meta.umap1, y=meta.umap2, c=c, s=1.3, rasterized=True)
plt.title(common.name_map_lpl_clusters["c7"], size=13)

despine(ax1)

# Plot the small umaps

genes = ["Eomes", "Cd27", "Il10", "Lag3", "Gzmk"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="grays", colors=["#dddddd", "000000"]
)

kwargs = dict(cmap=cmap, vmin=0, vmax=6, s=0.8, edgecolors="none")

for g, ax in zip(genes, small_axs):
    plt.sca(ax)

    exp = np.log2(expression.loc[g] + 1)
    exp = exp.loc[meta.index]

    pc = plt.scatter(x=meta.umap1, y=meta.umap2, c=exp, rasterized=True, **kwargs)

    despine(ax)
    plt.title(g, size=11)

plt.colorbar(
    pc, cax=cbar_ax, ticks=[0, 3, 6], label="Expression\n$log_2(CP10K+1)$"
)

plt.savefig('umap_grid.svg', dpi=300)


# %% Alternate - plot just a 3x2 grid:

fig = plt.figure(figsize=(8.5*1.2, 5*1.2))

gs = fig.add_gridspec(2, 3, left=0.1, right=0.85, top=0.9, bottom=0.1)

ax1 = fig.add_subplot(gs[0, 0])

ax1s = fig.add_subplot(gs[0, 1])
ax2s = fig.add_subplot(gs[0, 2])
ax3s = fig.add_subplot(gs[1, 0])
ax4s = fig.add_subplot(gs[1, 1])
ax5s = fig.add_subplot(gs[1, 2])

cbar_ax = fig.add_axes([0.88, 0.1, 0.01, 0.1])

small_axs = [ax1s, ax2s, ax3s, ax4s, ax5s]


def despine(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.sca(ax)
    plt.xticks([])
    plt.yticks([])


# Plot the big umap

plt.sca(ax1)

c_color = common.colors.cmap_lpl_clusters["c7"]

c = [c_color if x == "c7" else "#dddddd" for x in meta.Cluster]

plt.scatter(x=meta.umap1, y=meta.umap2, c=c, s=1.5, edgecolors="none", rasterized=True)
plt.title(common.name_map_lpl_clusters["c7"], size=13)

despine(ax1)

# Plot the small umaps

genes = ["Eomes", "Cd27", "Il10", "Lag3", "Gzmk"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    name="grays", colors=["#dddddd", "000000"]
)

kwargs = dict(cmap=cmap, vmin=0, vmax=6, s=1.5, edgecolors="none")

for g, ax in zip(genes, small_axs):
    plt.sca(ax)

    exp = np.log2(expression.loc[g] + 1)
    exp = exp.loc[meta.index]

    pc = plt.scatter(x=meta.umap1, y=meta.umap2, c=exp, **kwargs, rasterized=True)

    despine(ax)
    plt.title(g, size=11)

plt.colorbar(
    pc, cax=cbar_ax, ticks=[0, 3, 6], label="Expression\n$log_2(CP10K+1)$"
)

plt.savefig('umap_grid_alt.svg', dpi=300)
