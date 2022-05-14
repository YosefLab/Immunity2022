"""
Plot a larger, summary UMAP along with smaller, cluster-specific UMAPs
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import common

# %% Load Data

umap = pd.read_table("analysis/umap/umap.txt", index_col=0)
clusters = pd.read_table("analysis/cluster_together/clusters.txt", index_col=0)
meta = pd.read_table("analysis/data/meta.txt.gz", index_col=0)


data = meta.join(clusters).join(umap)
data = data.sample(frac=1)  # Randomize order for plots


def nobox():
    ax = plt.gca()
    for sp in ax.spines.values():
        sp.set_visible(False)

# %% Plot a Genotype UMAP

fig = plt.figure(figsize=(5, 4))

plt.scatter(
    x=data.umap1,
    y=data.umap2,
    s=1,
    color=[common.cmap_genotype[g] for g in data.Genotype],
    alpha=0.6,
    rasterized=True
)
plt.xticks([])
plt.yticks([])

# Create legend
plt.subplots_adjust(right=0.8)

handles = [
    Line2D(
        [], [], color=common.cmap_genotype['ctrl'],
        marker='o', markersize=5, linestyle=''),
    Line2D(
        [], [], color=common.cmap_genotype['KO'],
        marker='o', markersize=5, linestyle=''),
]

labels = [
    common.name_map_genotype['ctrl'], common.name_map_genotype['KO'],
]

plt.legend(handles, labels, loc='upper left', bbox_to_anchor=[.83, .95])
fig.patch.set_visible(False)
nobox()

# plt.show()
plt.savefig('umap_genotype.svg', dpi=300)

# %% Plot a Genotype UMAP again, but fully vectorized

fig = plt.figure(figsize=(5, 4))

plt.scatter(
    x=data.umap1,
    y=data.umap2,
    s=1,
    color=[common.cmap_genotype[g] for g in data.Genotype],
    alpha=0.6,
    rasterized=False
)
plt.xticks([])
plt.yticks([])

# Create legend
plt.subplots_adjust(right=0.8)

handles = [
    Line2D(
        [], [], color=common.cmap_genotype['ctrl'],
        marker='o', markersize=5, linestyle=''),
    Line2D(
        [], [], color=common.cmap_genotype['KO'],
        marker='o', markersize=5, linestyle=''),
]

labels = [
    common.name_map_genotype['ctrl'], common.name_map_genotype['KO'],
]

plt.legend(handles, labels, loc='upper left', bbox_to_anchor=[.83, .95])
fig.patch.set_visible(False)
nobox()

# plt.show()
plt.savefig('umap_genotype_vectorized.svg', dpi=300)
# %% Plot a single cluster - c2

cluster = "c2"
color = [
    common.cmap_lpl_clusters[x] if x == cluster else "#CCCCCC"
    for x in data.Cluster
]

fig = plt.figure(figsize=(2, 2))

plt.scatter(
    x=data.umap1,
    y=data.umap2,
    s=.3,
    color=color,
    alpha=0.6,
    rasterized=True
)
plt.xticks([])
plt.yticks([])

fig.patch.set_visible(False)
nobox()
plt.savefig('umap_cluster_{}.svg'.format(cluster), dpi=300)
# plt.show()

# %% Plot a single cluster - c9

cluster = "c9"
color = [
    common.cmap_lpl_clusters[x] if x == cluster else "#CCCCCC"
    for x in data.Cluster
]

fig = plt.figure(figsize=(2, 2))

plt.scatter(
    x=data.umap1,
    y=data.umap2,
    s=.3,
    color=color,
    alpha=0.6,
    rasterized=True
)
plt.xticks([])
plt.yticks([])

fig.patch.set_visible(False)
nobox()
plt.savefig('umap_cluster_{}.svg'.format(cluster), dpi=300)
# plt.show()

# %% Plot a single cluster - c7

cluster = "c7"
color = [
    common.cmap_lpl_clusters[x] if x == cluster else "#CCCCCC"
    for x in data.Cluster
]

fig = plt.figure(figsize=(2, 2))

plt.scatter(
    x=data.umap1,
    y=data.umap2,
    s=.3,
    color=color,
    alpha=0.6,
    rasterized=True
)
plt.xticks([])
plt.yticks([])

fig.patch.set_visible(False)
nobox()
plt.savefig('umap_cluster_{}.svg'.format(cluster), dpi=300)
# plt.show()

# %% Plot All Clusters

color = [
    common.cmap_lpl_clusters[x]
    for x in data.Cluster
]

fig = plt.figure(figsize=(5.5, 4))

plt.scatter(
    x=data.umap1,
    y=data.umap2,
    s=.6,
    color=color,
    alpha=0.6,
    rasterized=True
)
plt.xticks([])
plt.yticks([])

fig.patch.set_visible(False)
nobox()

# Create legend
handles = [
    Line2D(
        [], [], color=common.cmap_lpl_clusters[x],
        marker='o', markersize=5, linestyle='')
    for x in common.order_lpl_clusters
]

labels = [
    common.name_map_lpl_clusters[x] for x in common.order_lpl_clusters
]

plt.subplots_adjust(right=0.7)
fig.legend(handles, labels, bbox_to_anchor=[0.7, .9], loc='upper left')

plt.savefig('umap_cluster_all.svg'.format(cluster), dpi=300)
