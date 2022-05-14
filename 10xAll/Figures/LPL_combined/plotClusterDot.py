import feather
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import bio_utils.plots
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colors import to_rgb, rgb_to_hsv
import common
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

# %% Import data


clusters = pd.read_table("analysis/cluster_together/clusters.txt", index_col=0)

scaled = feather.read_dataframe("analysis/data/expression_scaled.feather") \
    .set_index('index')

# %% Cluster percent expression and cluster mean

cluster_means = scaled.T.join(clusters).groupby("Cluster").mean() / 3000 * 10000
cluster_means.columns.name = 'Gene'
cluster_means = cluster_means.T
cluster_means_norm = cluster_means \
    .subtract(cluster_means.mean(axis=1), axis=0) \
    .divide(cluster_means.std(axis=1), axis=0) \


cluster_percents = (scaled > 0).T.join(clusters).groupby("Cluster").mean() * 100
cluster_percents.columns.name = 'Gene'
cluster_percents = cluster_percents.T

# %% plot function

dot_genes_c2 = [
    "Ifitm1",
    "Itgb1",
    "Ly6c2",
    "S100a4",
    "Lgals3",
    "Klrd1",
    "Klf2",
    "Pglyrp1",
    "Anxa1",
    "S100a6",
    "Crip1",
    "Gzma",
    "Lgals1",
    "Ifitm3",
    "Klrk1",
    "AA467197",
    "Emp3",
    "Ifitm2",
    "S1pr1",
    "Vim",
]

dot_genes_c9 = [
    "Cd160",
    "Cd7",
    "Itgae",
    "Gzma",
    "Ccr9",
    "1810041H14Rik",
    "Sash3",
    "Cd38",
    "Abi3",
    "Sh2d3c",
    "Tbc1d10c",
    "AW112010",
    "Cd226",
    "Ptpn7",
    "Cd48",
    "Susd3",
    "Gzmk",
    "Dtx1",
    "Rinl",
    "Rsrp1",
]


def make_cluster_dot(
    data, x, y,
    size,
    hue,
    cmap, vmin, vmax,
    legend_sizes,
    x_order=None,
    y_order=None,
    size_range=None,
    bbox_cax=[.65, .55, .02, .2], bbox_legend=(.6, .5)
):

    size_range = None
    if size_range is None:
        size_range = (data[size].min(), data[size].max())

    if x_order is None:
        x_order = np.unique(data[x])
    if y_order is None:
        y_order = np.unique(data[y])

    ax, x_map, y_map, max_circle_area = bio_utils.plots.create_dot_plot_axes(
        data=data,
        x=x, x_order=x_order,
        y=y, y_order=y_order,
        bpad=0.2, lpad=0.2, square_pad=0.8
    )

    def size_scale_fn(x):
        scale = max_circle_area / (size_range[1] - size_range[0])
        offset = size_range[0]
        return scale*(x-offset)

    if isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap)
    else:
        cmap_obj = cmap

    norm = Normalize(vmin, vmax, clip=True)

    dd = data.copy()
    dd[x] = dd[x].map(x_map)
    dd[y] = dd[y].map(y_map)
    dd[size] = [size_scale_fn(x) for x in dd[size]]
    dd[hue] = [cmap_obj(norm(x)) for x in dd[hue]]

    plt.scatter(
        x=dd[x],
        y=dd[y],
        c=dd[hue],
        s=dd[size],
        linewidths=.5,
        edgecolors='black',
    )

    bio_utils.plots.add_size_legend(
        sizes=legend_sizes, size_to_pts_fn=size_scale_fn,
        bbox_to_anchor=bbox_legend, labelspacing=.4, title='% Detected',
        edgecolor='white', loc='upper left', fancybox=False, mew=.5
    )

    fig = plt.gcf()
    cax = fig.add_axes(bbox_cax)
    fig.colorbar(ScalarMappable(norm, cmap_obj), cax=cax, label='logFC')

    ax.set_xticklabels(
        [common.name_map_lpl_clusters[x.get_text()] for x in ax.get_xticklabels()]
    )
    plt.sca(ax)
    plt.xticks(rotation=90)


# %%  Plot c9

plt.figure()

dot_genes = dot_genes_c9

cmeans = cluster_means_norm.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Mean")

cper = cluster_percents.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Percent")

data = cmeans.merge(cper)

make_cluster_dot(
    data=data,
    x="Cluster",
    y="Gene",
    size="Percent",
    hue="Mean",
    cmap="RdYlBu_r",
    vmin=-2,
    vmax=2,
    legend_sizes=[25, 50, 75, 100],
    x_order=common.order_lpl_clusters,
    y_order=dot_genes[::-1],
    size_range=None,
    bbox_cax=[.65, .55, .02, .2],
    bbox_legend=(.6, .5),
)

# plt.show()
plt.savefig('c9_dot.svg')

# %%  Plot c2

plt.figure()


dot_genes = dot_genes_c2


cmeans = cluster_means_norm.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Mean")

cper = cluster_percents.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Percent")

data = cmeans.merge(cper)

make_cluster_dot(
    data=data,
    x="Cluster",
    y="Gene",
    size="Percent",
    hue="Mean",
    cmap="RdYlBu_r",
    vmin=-2,
    vmax=2,
    legend_sizes=[25, 50, 75, 100],
    x_order=common.order_lpl_clusters,
    y_order=dot_genes[::-1],
    size_range=None,
    bbox_cax=[.65, .55, .02, .2],
    bbox_legend=(.6, .5),
)

# plt.show()
plt.savefig('c2_dot.svg')

# %%

##############################
# Alternate Color Scale      #
##############################

# Here, each row is linearly scaled with the top gene
# representing the maximum value (per cluster)


def make_cluster_dot_v2(
    data, x, y,
    size,
    hue,
    cmap, vmin, vmax,
    legend_sizes,
    x_order=None,
    y_order=None,
    size_range=None,
    bbox_cax=[.65, .55, .02, .2], bbox_legend=(.6, .5)
):

    size_range = None
    if size_range is None:
        size_range = (data[size].min(), data[size].max())

    if x_order is None:
        x_order = np.unique(data[x])
    if y_order is None:
        y_order = np.unique(data[y])

    ax, x_map, y_map, max_circle_area = bio_utils.plots.create_dot_plot_axes(
        data=data,
        x=x, x_order=x_order,
        y=y, y_order=y_order,
        bpad=0.2, lpad=0.2, square_pad=0.8
    )

    def size_scale_fn(x):
        scale = max_circle_area / (size_range[1] - size_range[0])
        offset = size_range[0]
        return scale*(x-offset)

    if isinstance(cmap, str):
        cmap_obj = plt.get_cmap(cmap)
    else:
        cmap_obj = cmap

    norm = Normalize(vmin, vmax, clip=True)

    dd = data.copy()
    dd[x] = dd[x].map(x_map)
    dd[y] = dd[y].map(y_map)
    dd[size] = [size_scale_fn(x) for x in dd[size]]
    dd[hue] = [cmap_obj(norm(x)) for x in dd[hue]]

    plt.scatter(
        x=dd[x],
        y=dd[y],
        c=dd[hue],
        s=dd[size],
        linewidths=.5,
        edgecolors='black',
    )

    bio_utils.plots.add_size_legend(
        sizes=legend_sizes, size_to_pts_fn=size_scale_fn,
        bbox_to_anchor=bbox_legend, labelspacing=.4, title='% Detected',
        edgecolor='white', loc='upper left', fancybox=False, mew=.5
    )

    fig = plt.gcf()
    cax = fig.add_axes(bbox_cax)
    fig.colorbar(ScalarMappable(norm, cmap_obj), cax=cax, label='Relative\nExpression', ticks=[0, 1])

    ax.set_xticklabels(
        [common.name_map_lpl_clusters[x.get_text()] for x in ax.get_xticklabels()]
    )
    plt.sca(ax)
    plt.xticks(rotation=90)

# %%  Plot c9

plt.figure()

dot_genes = dot_genes_c9

cmaxes = cluster_means.loc[dot_genes]
cmaxes = cmaxes.divide(cmaxes.max(axis=1), axis=0) \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Max")

cper = cluster_percents.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Percent")

data = cmaxes.merge(cper)

hue = rgb_to_hsv(to_rgb(common.cmap_lpl_clusters['c9']))[0]*3 + .85
cmap = sns.cubehelix_palette(start=hue, light=1, rot=0, as_cmap=True)

make_cluster_dot_v2(
    data=data,
    x="Cluster",
    y="Gene",
    size="Percent",
    hue="Max",
    cmap=cmap,
    vmin=0,
    vmax=1,
    legend_sizes=[25, 50, 75, 100],
    x_order=common.order_lpl_clusters,
    y_order=dot_genes[::-1],
    size_range=None,
    bbox_cax=[.65, .55, .02, .2],
    bbox_legend=(.6, .5),
)

plt.show()

# %%  Plot c2

plt.figure()

dot_genes = dot_genes_c2

cmaxes = cluster_means.loc[dot_genes]
cmaxes = cmaxes.divide(cmaxes.max(axis=1), axis=0) \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Max")

cper = cluster_percents.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Percent")

data = cmaxes.merge(cper)

hue = rgb_to_hsv(to_rgb(common.cmap_lpl_clusters['c2']))[0]*3 + .85
cmap = sns.cubehelix_palette(start=hue, light=1, rot=0, as_cmap=True)

make_cluster_dot_v2(
    data=data,
    x="Cluster",
    y="Gene",
    size="Percent",
    hue="Max",
    cmap=cmap,
    vmin=0,
    vmax=1,
    legend_sizes=[25, 50, 75, 100],
    x_order=common.order_lpl_clusters,
    y_order=dot_genes[::-1],
    size_range=None,
    bbox_cax=[.65, .55, .02, .2],
    bbox_legend=(.6, .5),
)

plt.show()

##############################
# Plot T-cell Subset Genes   #
##############################

# %%

plt.figure()

dot_genes = [
    "Tbx21",
    "Ifng",
    "Il17a",
    "Il17f",
    "Il21r",
    "Il2rg",
    "Il12rb1",
    "Il12rb2",
    "Il23r",
    "Stat3",
    "Rorc",
    "Stat1",
    "Stat4",
    "Ccr6",
    "Tnf",
    "Csf2",
    "Il22",
    "Foxp3",
    "Klrb1",
    "Prdm1",
    "Il2ra",
    "Il2rb",
    "Ccr4",
    "Cxcr3",
]


cmeans = cluster_means_norm.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Mean")

cper = cluster_percents.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Percent")

data = cmeans.merge(cper)

make_cluster_dot(
    data=data,
    x="Cluster",
    y="Gene",
    size="Percent",
    hue="Mean",
    cmap="RdYlBu_r",
    vmin=-2,
    vmax=2,
    legend_sizes=[25, 50, 75, 100],
    x_order=common.order_lpl_clusters,
    y_order=dot_genes[::-1],
    size_range=None,
    bbox_cax=[.65, .55, .02, .2],
    bbox_legend=(.6, .5),
)

plt.show()

# %% Try this with c7 too

dot_genes_c7 = [
    'Eomes',
    'Gzmk',
    'Cd27',
    'Lag3',
    'Pdcd1',
    'Cxcr5',
    'Il10',
    'Il10ra',
    'Bend4',
    'Esm1',
    'Chchd10',
    'Klrg1',
    'Gpr18',
    'Ccr5',
    'Irf7',
    'Sema4a',
    'Rgs16',
    'Klf3',
    'Sh2d1a',
    'Nr4a2',
]

plt.figure()


dot_genes = dot_genes_c7

cmeans = cluster_means_norm.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Mean")

cper = cluster_percents.loc[dot_genes] \
    .reset_index() \
    .melt(id_vars="Gene", value_name="Percent")

data = cmeans.merge(cper)

make_cluster_dot(
    data=data,
    x="Cluster",
    y="Gene",
    size="Percent",
    hue="Mean",
    cmap="RdYlBu_r",
    vmin=-2,
    vmax=2,
    legend_sizes=[25, 50, 75, 100],
    x_order=common.order_lpl_clusters,
    y_order=dot_genes[::-1],
    size_range=None,
    bbox_cax=[.65, .55, .02, .2],
    bbox_legend=(.6, .5),
)

# plt.show()
plt.savefig('c7_dot.svg')
