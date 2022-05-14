import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bio_utils.plots
import seaborn as sns
import feather
import common

plt.rcParams['svg.fonttype'] = 'none'

# %% Load counts to get averages...

scaled = feather.read_dataframe("analysis/data/expression_scaled.feather").set_index('index')
meta = pd.read_table("analysis/data/meta.txt.gz", index_col=0)
gene_info = pd.read_table("analysis/data/genes.tsv", index_col=0)

tissue_means = scaled.T.join(meta.Tissue).groupby('Tissue').mean().T / 3000 * 10000
tissue_means.index.name = 'GeneSymbol'
gene_means = scaled.mean(axis=1) / 3000 * 10000

# %% Load DE results for each subset
tissue_dir = {
    'SI_LPL': 'combined_SI_LPL_main',
    'Spleen': 'combined_Spleen_main',
    'Colon_LPL': 'combined_Colon_LPL_main',
    'Colon_IEL': 'combined_Colon_IEL_main',
}


geno_de = {}
for tiss in tissue_dir:
    geno_de_tiss = pd.read_table(
        "../../{}/de_all/de_results_edger.txt.gz".format(tissue_dir[tiss]),
    )
    geno_de_tiss['Tissue'] = tiss
    geno_de[tiss] = geno_de_tiss

geno_de_all = pd.concat((geno_de.values()), axis=0)


# Number of DE genes per tissue?
all_geno_genes = set()
for x in geno_de:
    selected = (geno_de[x].FDR < .05) & (geno_de[x].logFC.abs() > .5)
    genes = geno_de[x][selected].GeneSymbol
    print(x, len(genes))
    all_geno_genes |= set(genes)

# %% Make a heatmap

data_x = geno_de_all.loc[
    [x in all_geno_genes for x in geno_de_all.GeneSymbol]
]

data = data_x.pivot(index='GeneSymbol', columns='Tissue', values='logFC')
data = data.fillna(0)  # Genes are NA in Tissue if they were pre-filtered

data_fdr = data_x.pivot(index='GeneSymbol', columns='Tissue', values='FDR')
data_fdr = data_fdr.fillna(1)  # Genes are NA in Tissue if they were pre-filtered

# %% Dot plot of specific genes

dot_genes1 = [
    # Up everywhere in KO
    'Dnajb1',
    'Hspa1a',
    'Cxcl10',
    'Cxcr5',
    # Upregulated in ctrl everywhere (generally)
    'S1pr1',
    'Gzma',
    'Dusp2',
    'Gpr18',
    'Lgals3',
    'Cd160',
    'Zfp36l2',
    'Il7r',
    'Itgb7',
    'Itgae',
    'Ccr9',
]

dot_genes2 = [
    # Upregulated in ctrl in periphery
    'S100a4',
    'Ccl5',
    # Up in LPL in KO
    'Ahr',
    'Socs3',
    # Up everywhere except spleen in KO
    'Ccl20',
    'Cd69',
    'Ctla4',
    'Fos',
    'Jun',
    'Fosb',
]


logfcs = (data*-1).reset_index().melt(id_vars='GeneSymbol', value_name='logFC')
sig = (data_fdr < .05).reset_index().melt(id_vars='GeneSymbol', value_name='Significant')
sizes = np.log2(tissue_means+1).reset_index().melt(id_vars='GeneSymbol', value_name='Expression')

dot_data = logfcs.merge(sig, on=["GeneSymbol", "Tissue"]) \
    .merge(sizes, on=["GeneSymbol", "Tissue"])


# %%
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def make_dot_plot(data, x, y, size,
                  hue, cmap, vmin, vmax,
                  legend_sizes,
                  x_order=None, y_order=None,
                  size_range=None,
                  bbox_cax=[.41, .55, .02, .2], bbox_legend=(.4, .5)):

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
        bpad=0.2, lpad=0.2
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

    dd_x = dd.loc[dd.Significant==False]

    plt.scatter(
        x=dd_x[x],
        y=dd_x[y],
        c=dd_x[hue],
        s=dd_x[size],
    )


    dd_x = dd.loc[dd.Significant==True]

    plt.scatter(
        x=dd_x[x],
        y=dd_x[y],
        c=dd_x[hue],
        s=dd_x[size],
        linewidths=1,
        edgecolors='black'
    )

    bio_utils.plots.add_size_legend(
        sizes=legend_sizes, size_to_pts_fn=size_scale_fn,
        bbox_to_anchor=bbox_legend, labelspacing=1.2, title='Expression',
        edgecolor='white', loc='upper left', fancybox=False
    )

    fig = plt.gcf()
    cax = fig.add_axes(bbox_cax)
    fig.colorbar(ScalarMappable(norm, cmap_obj), cax=cax, label='logFC')

    ax.set_xticklabels(
        [common.name_map_tissue[x.get_text()] for x in ax.get_xticklabels()]
    )
    plt.sca(ax)
    plt.xticks(rotation=90)

# %%

from matplotlib.colors import to_rgb, rgb_to_hsv

cmap = "RdBu_r"
# cmap = sns.diverging_palette(
#     h_neg=rgb_to_hsv(to_rgb(common.cmap_genotype['KO']))[0]*360,
#     h_pos=rgb_to_hsv(to_rgb(common.cmap_genotype['ctrl']))[0]*360,
#     s=60, l=60, as_cmap=True
# )

dot_data_sub = dot_data.loc[
    [x in dot_genes1 for x in dot_data.GeneSymbol]
]

fig = plt.figure(figsize=(5, 7))
make_dot_plot(
    dot_data_sub,
    x="Tissue",
    y="GeneSymbol",
    size="Expression",
    hue="logFC",
    cmap=cmap,
    vmin=-1,
    vmax=1,
    legend_sizes=[1, 3, 5],
    x_order=common.order_tissue,
    y_order=dot_genes1,
    bbox_cax=[.56, .55, .02, .2],
    bbox_legend=(.55, .5)
)
# plt.show()
plt.savefig("tissue_geno_dot1.svg")

# %%

dot_data_sub = dot_data.loc[
    [x in dot_genes2 for x in dot_data.GeneSymbol]
]

fig = plt.figure(figsize=(5, 5))
make_dot_plot(
    dot_data_sub,
    x="Tissue",
    y="GeneSymbol",
    size="Expression",
    hue="logFC",
    cmap=cmap,
    vmin=-1,
    vmax=1,
    legend_sizes=[1, 3, 5],
    x_order=common.order_tissue,
    y_order=dot_genes2,
    bbox_cax=[.56, .55, .02, .2],
    bbox_legend=(.55, .5)
)
#plt.show()
plt.savefig("tissue_geno_dot2.svg")
