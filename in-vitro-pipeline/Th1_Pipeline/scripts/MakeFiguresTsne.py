# %%
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib import gridspec, cm
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

#tpm = pd.read_table("normalized_noqpc/tpm.txt.gz", index_col=0)
tpm = pd.read_table("data_filtered/tpm_all.txt.gz", index_col=0)
meta = pd.read_table("data_filtered/meta.txt.gz", index_col=0)
tsne = pd.read_table("tsne/tsne.txt", index_col=0)
clusters = pd.read_table("clusters/clusters.txt", index_col=0)


# %% Meta-data Plot

fig = plt.figure(figsize=(4, 4))

groups = [
    ('ctrl', 'eGFP_neg', 'black', 'o'),
    ('ctrl', 'eGFP_pos', '#287b2e', 'o'),
    ('KO', 'eGFP_neg', 'black', 'x'),
    ('KO', 'eGFP_pos', '#287b2e', 'x'),
]

label_map = {
    'ctrl': '$IL23R^{WT/eGFP}$',
    'KO': '$IL23R^{eGFP/eGFP}$',
    'eGFP_neg': 'GFP-',
    'eGFP_pos': 'GFP+',
}


for geno, gfp, color, marker in groups:

    ii = meta.index[
        (meta.Genotype == geno) &
        (meta.GFP == gfp)
    ]

    plt.plot(tsne.tsne1[ii], tsne.tsne2[ii], ms=3.5,
             color=color, marker=marker, linestyle='None',
             label="{} {}".format(label_map[gfp], label_map[geno])
             )

plt.xlabel('TSNE-1')
plt.ylabel('TSNE-2')
plt.xticks([])
plt.yticks([])
plt.legend()
fig.patch.set_visible(False)
plt.savefig('Figures/tsne_meta.svg')

# %% Plot Gene1

fig = plt.figure(figsize=(3.5, 3))

gene = 'IL2'
cmap = LinearSegmentedColormap.from_list("grays", ["#DDDDDD", "#222222"])

x = tsne.tsne1
y = tsne.tsne2
c = np.log1p(tpm.loc[gene][tsne.index])

plt.scatter(x, y, c=c, s=3, cmap=cmap)
plt.title(gene)

plt.xticks([])
plt.yticks([])
plt.colorbar(shrink=0.5)
fig.patch.set_visible(False)
plt.savefig('Figures/tsne_{}.svg'.format(gene))
# plt.show()

# %%

fig = plt.figure(figsize=(3.5, 3))

gene = 'EGR2'
cmap = LinearSegmentedColormap.from_list("grays", ["#CCCCCC", "#000000"])

x = tsne.tsne1
y = tsne.tsne2
c = np.log1p(tpm.loc[gene][tsne.index])

plt.scatter(x, y, c=c, s=3, cmap=cmap)
plt.title(gene)

plt.xticks([])
plt.yticks([])
plt.colorbar(shrink=0.5)
fig.patch.set_visible(False)
plt.savefig('Figures/tsne_{}.svg'.format(gene))
#plt.show()

# %%

fig = plt.figure(figsize=(3.5, 3))

gene = 'NR4A1'
cmap = LinearSegmentedColormap.from_list("grays", ["#CCCCCC", "#000000"])

x = tsne.tsne1
y = tsne.tsne2
c = np.log1p(tpm.loc[gene][tsne.index])

plt.scatter(x, y, c=c, s=3, cmap=cmap)
plt.title(gene)

plt.xticks([])
plt.yticks([])
plt.colorbar(shrink=0.5)
fig.patch.set_visible(False)
plt.savefig('Figures/tsne_{}.svg'.format(gene))
#plt.show()

# %% Plot clusters

fig = plt.figure(figsize=(2.85, 3))

x = tsne.tsne1
y = tsne.tsne2
colors = plt.get_cmap('tab20').colors

clusters = clusters.loc[tsne.index]

ii = clusters.Cluster == 'c5'
plt.plot(x[~ii], y[~ii], 'o', ms=3**.5, color="#CCCCCC")
plt.plot(x[ii], y[ii], 'o', ms=3**.5, color=colors[0], label="Cluster 5")
plt.title("Clusters")

plt.xticks([])
plt.yticks([])
plt.legend()
fig.patch.set_visible(False)
plt.savefig('Figures/tsne_c5.svg')

# plt.show()

# %% Cluster-proportion plot

cluster_id = 'c0'

plot_data = meta.copy().join(clusters)
plot_data['x'] = 0
plot_data.loc[plot_data.Cluster == cluster_id, 'x'] = 1

cluster_counts = plot_data.groupby(['Genotype', 'GFP'])['x'].sum()
cluster_denom = plot_data.groupby(['Genotype', 'GFP'])['x'].size()

cluster_proportion = cluster_counts / cluster_denom * 100
cluster_proportion = cluster_proportion.reset_index()

groups = [
    ('ctrl', 'eGFP_neg'),
    ('ctrl', 'eGFP_pos'),
    ('KO', 'eGFP_neg'),
    ('KO', 'eGFP_pos'),
]

label_map = {
    'ctrl': '$IL23R^{WT/eGFP}$',
    'KO': '$IL23R^{eGFP/eGFP}$',
    'eGFP_neg': 'GFP-',
    'eGFP_pos': 'GFP+',
}

x = []
height = []

for geno, gfp in groups:

    x.append("{} {}".format(label_map[gfp], label_map[geno]))
    height.append(cluster_proportion.loc[
        (cluster_proportion.Genotype == geno) &
        (cluster_proportion.GFP == gfp),
        'x'].values[0])

fig = plt.figure(figsize=(3, 3))
plt.grid(axis='y')
plt.gca().set_axisbelow(True)
plt.bar(x=x, height=height)
plt.xlabel('Sample')
plt.ylabel('Proportion (%)')
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.title('{} Composition'.format(cluster_id))
fig.patch.set_visible(False)
plt.savefig('Figures/cluster_{}_proportion.svg'.format(cluster_id))

# %% Cluster-proportion plot

cluster_id = 'c5'

plot_data = meta.copy().join(clusters)
plot_data['x'] = 0
plot_data.loc[plot_data.Cluster == cluster_id, 'x'] = 1

cluster_counts = plot_data.groupby(['Genotype', 'GFP'])['x'].sum()
cluster_denom = plot_data.groupby(['Genotype', 'GFP'])['x'].size()

cluster_proportion = cluster_counts / cluster_denom * 100
cluster_proportion = cluster_proportion.reset_index()

groups = [
    ('ctrl', 'eGFP_neg'),
    ('ctrl', 'eGFP_pos'),
    ('KO', 'eGFP_neg'),
    ('KO', 'eGFP_pos'),
]

label_map = {
    'ctrl': '$IL23R^{WT/eGFP}$',
    'KO': '$IL23R^{eGFP/eGFP}$',
    'eGFP_neg': 'GFP-',
    'eGFP_pos': 'GFP+',
}

x = []
height = []

for geno, gfp in groups:

    x.append("{} {}".format(label_map[gfp], label_map[geno]))
    height.append(cluster_proportion.loc[
        (cluster_proportion.Genotype == geno) &
        (cluster_proportion.GFP == gfp),
        'x'].values[0])

fig = plt.figure(figsize=(3, 3))
plt.grid(axis='y')
plt.gca().set_axisbelow(True)
plt.bar(x=x, height=height)
plt.xlabel('Sample')
plt.ylabel('Proportion (%)')
plt.subplots_adjust(left=0.2, bottom=0.2)
plt.title('{} Composition'.format(cluster_id))
fig.patch.set_visible(False)
plt.savefig('Figures/cluster_{}_proportion.svg'.format(cluster_id))
#plt.show()

# %% TSNE Grid

genes = ['IFNG', 'IL23R', 'TBX21',
         'IL17F', 'RORC', 'CCR6',
         'TNF', 'IL2', 'CCL4']

fig, axs = plt.subplots(3, 3, figsize=(7,7))

cmap = LinearSegmentedColormap.from_list("grays", ["#DDDDDD", "#222222"])

for g, ax in zip(genes, axs.ravel()):

    plt.sca(ax)

    for sp in ax.spines:
        ax.spines[sp].set_visible(False)

    x = tsne.tsne1
    y = tsne.tsne2
    c = np.log1p(tpm.loc[g][tsne.index])

    sc = plt.scatter(x, y, c=c, s=3,
                cmap=cmap, vmin=0, vmax=6,
                edgecolors=None, rasterized=True)
    plt.title(g)

    plt.xticks([])
    plt.yticks([])

plt.subplots_adjust(hspace=0.3)

# Add a colorbar
cax = fig.add_axes(
    [.94, .02, .01, .15],
    label='cax'
)
plt.colorbar(sc, cax=cax)
plt.savefig('Figures/tsne_grid.svg', dpi=300)

# plt.show()

# %% TSNE Grid - more genes

genes = ['IFNG', 'IL12RB2', 'TBX21', 'CXCR3',
         'RORC', 'IL17A', 'IL17F', 'CCR6',
         'TNF', 'IL2', 'IL23R', 'IL21R']

fig, axs = plt.subplots(3, 4, figsize=(8,7))

cmap = LinearSegmentedColormap.from_list("grays", ["#DDDDDD", "#222222"])

for g, ax in zip(genes, axs.ravel()):

    plt.sca(ax)

    if g == '':
        ax.remove()
        continue

    for sp in ax.spines:
        ax.spines[sp].set_visible(False)

    x = tsne.tsne1
    y = tsne.tsne2
    c = np.log1p(tpm.loc[g][tsne.index])

    sc = plt.scatter(x, y, c=c, s=3,
                cmap=cmap, vmin=0, vmax=6,
                edgecolors=None, rasterized=True)
    plt.title(g)

    plt.xticks([])
    plt.yticks([])

plt.subplots_adjust(hspace=0.3)

# Add a colorbar
cax = fig.add_axes(
    [.96, .11, .01, .15],
    label='cax'
)
plt.colorbar(sc, cax=cax)
plt.sca(cax)
plt.ylabel('Log(TPM+1)', labelpad=-35)
plt.savefig('Figures/tsne_grid2.svg', dpi=300)

# plt.show()

# %% Instead of tSNE grid, do a dot plot

import bio_utils
import matplotlib.colors
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def make_dot(
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

    if size_range is None:
        size_range = (data[size].min(), data[size].max())

    if x_order is None:
        x_order = np.unique(data[x])
    if y_order is None:
        y_order = np.unique(data[y])
    y_order = y_order[::-1]

    print(x_order)
    print(y_order)

    ax, x_map, y_map, max_circle_area = bio_utils.plots.create_dot_plot_axes(
        data=data,
        x=x, x_order=x_order,
        y=y, y_order=y_order,
        bpad=0.2, lpad=0.25, square_pad=0.8
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
        bbox_to_anchor=bbox_legend, labelspacing=1, title='Detection (%)',
        edgecolor='white', loc='lower left', fancybox=False, mew=.5
    )

    fig = plt.gcf()
    cax = fig.add_axes(bbox_cax)
    fig.colorbar(ScalarMappable(norm, cmap_obj), cax=cax,
                 orientation='horizontal', ticks=[0, vmax/2, vmax]
                 )
    # cax.set_yticklabels(['0', 'Max'])
    cax.set_xlabel('Mean Expression\n$log(TPM+1)$', labelpad=15)
    cax.xaxis.set_label_position('top')

    # ax.set_xticklabels(
    #     [common.name_map_lpl_clusters[x.get_text()] for x in ax.get_xticklabels()]
    # )
    plt.sca(ax)
    plt.xticks(rotation=90)


genes = ['IFNG', 'IL12RB2', 'TBX21', 'CXCR3',
         'RORC', 'IL17A', 'IL17F', 'CCR6',
         'TNF', 'IL2', 'IL23R', 'IL21R']

c = np.log1p(tpm.loc[genes]).T
c = c.join(meta[['Genotype', 'GFP']])
c.index.name = 'Gene'
c = c.groupby(['Genotype', 'GFP']).mean().reset_index()
c = c.melt(id_vars=['Genotype', 'GFP'], var_name='Gene', value_name='Mean Expression')

d = (tpm.loc[genes] > 0).T
d = d.join(meta[['Genotype', 'GFP']])
d.index.name = 'Gene'
d = (d.groupby(['Genotype', 'GFP']).mean() * 100).reset_index()
d = d.melt(id_vars=['Genotype', 'GFP'], var_name='Gene', value_name='Detection')

plot_data = c.merge(d, on=['Genotype', 'GFP', 'Gene'])
plot_data['Sample'] = ['{}-{}'.format(x, y) for x, y in zip(plot_data.Genotype, plot_data.GFP)]

sample_map = {
    'KO-eGFP_neg': '$IL23R^{eGFP/eGFP}$, $GFP-$',
    'ctrl-eGFP_neg': '$IL23R^{wt/eGFP}$, $GFP-$',
    'KO-eGFP_pos': '$IL23R^{eGFP/eGFP}$, $GFP+$',
    'ctrl-eGFP_pos': '$IL23R^{wt/eGFP}$, $GFP+$',
}
plot_data['Sample'] = plot_data['Sample'].map(sample_map)

legend_sizes = [25, 50, 75]
cmap = sns.cubehelix_palette(start=2.63, light=1, rot=0, as_cmap=True, hue=1)
# cmap = "RdBu_r"
vmin = 0
vmax = 6

plt.figure(figsize=(7, 5))

y_order = [
    '$IL23R^{eGFP/eGFP}$, $GFP-$',
    '$IL23R^{wt/eGFP}$, $GFP-$',
    '$IL23R^{eGFP/eGFP}$, $GFP+$',
    '$IL23R^{wt/eGFP}$, $GFP+$',
]

make_dot(
    data=plot_data,
    x="Gene",
    y="Sample",
    size="Detection",
    hue="Mean Expression",
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
    legend_sizes=legend_sizes,
    x_order=genes,
    y_order=y_order,
    size_range=None,
    bbox_cax=[0.5, 0.65, 0.15, 0.02],
    bbox_legend=(0.2, 0.6),
)

ax = plt.gca()
ax.tick_params(axis='y', left=False)

plt.savefig('Figures/marker_dots.svg')
# plt.show()
