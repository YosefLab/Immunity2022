import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list

plt.rcParams['svg.fonttype'] = 'none'

# %% Load the DE Data


data = pd.read_table("analysis/cluster_together/cluster_de_genotype_edgeR/cluster_de_results.txt.gz")


# gwas = pd.read_table("markers_from_mathias.txt").iloc[:, 0].tolist()
# gwas = [x.lower() for x in gwas]

gwas = open("../../../GWAS/genes_delange_implicated.txt").readlines()
gwas = [x.strip().lower() for x in gwas]


data_sub = data.loc[
    (data.cluster.isin({'c7', 'c9', 'c2'})) &
    (data.FDR < .1) &
    (data.GeneSymbol.str.lower().isin(gwas))
]


genes_to_plot = data_sub.GeneSymbol.unique()

# Hmm, have to figure out the order
# %% Plot a heatmap


plot_data = data.loc[
    (data.cluster.isin({'c7', 'c9', 'c2'})) &
    (data.GeneSymbol.isin(genes_to_plot))
]

annot_data = pd.pivot(
    plot_data, index='GeneSymbol',
    columns='cluster', values='FDR'
)

annot_data = annot_data < .1
annot_data[annot_data] = '*'
annot_data[annot_data != '*'] = ''

plot_data = pd.pivot(
    plot_data, index='GeneSymbol',
    columns='cluster', values='logFC'
) * -1

plot_data = plot_data[['c2', 'c9', 'c7']]
annot_data = annot_data.loc[plot_data.index].loc[:, plot_data.columns]

# %% Not bad.  Still could do better with the ordering

x = pd.crosstab(data_sub.GeneSymbol, data_sub.cluster)

g1 = x.index[
    (x.c2 == 1) & (x.c7 == 0) & (x.c9 == 0)
]

g2 = x.index[
    (x.c2 == 0) & (x.c7 == 0) & (x.c9 == 1)
]

g3 = x.index[
    (x.c2 == 0) & (x.c7 == 1) & (x.c9 == 0)
]

g4 = x.index[
    (x.c2 == 1) & (x.c7 == 1) & (x.c9 == 0)
]

g5 = x.index[
    (x.c2 == 0) & (x.c7 == 1) & (x.c9 == 1)
]

g6 = x.index[
    (x.c2 == 1) & (x.c7 == 0) & (x.c9 == 1)
]

g7 = x.index[
    (x.c2 == 1) & (x.c7 == 1) & (x.c9 == 1)
]

groups = [g1, g2, g3, g4, g5, g6, g7]


def clust_group(genes):

    if len(genes) < 2:
        return genes

    dd = plot_data.loc[genes]
    ii = leaves_list(
        linkage(dd, metric='euclidean', method='average')
    )

    return dd.index[ii]

groups = [clust_group(g) for g in groups]

gene_plot_order = []
for group in groups:
    gene_plot_order.extend(group.tolist())

plot_data = plot_data.loc[gene_plot_order]
annot_data = annot_data.loc[plot_data.index].loc[:, plot_data.columns]

fig = plt.figure(figsize=(5, 10))

ax = fig.add_axes([.2, .1, .5, .8])
cb_ax = fig.add_axes([.8, .1, .02, .1])

plt.sca(ax)
hm = sns.heatmap(
    plot_data, cmap='RdBu_r', yticklabels=True, annot=annot_data, fmt='s',
    vmin=-1, vmax=1,
    ax=ax, cbar_ax=cb_ax, cbar_kws=dict(
        ticks=[-1, 0, 1],
    ),
)
plt.xlabel('')
plt.ylabel('')
plt.sca(cb_ax)
plt.ylabel('Ctrl vs. KO\n$log_2$ Fold-Change')
# plt.show()
plt.savefig('gwas_heatmap_c2c7c9.svg')


# %% How do these intersect with the Smillie data?
# What are the GWAS genes DE in their CD4 subset?

# regev_de = pd.read_excel(
#     "../../../../externalData/Smillie2019/supplements/1-s2.0-S0092867419307329-mmc4.xlsx",
#     sheet_name="Adaptive (Inflamed vs. Healthy)"
# )
# 
# all_res = []
# for ss in regev_de.ident.unique():
#     cd4_de = regev_de.loc[regev_de.ident == ss]
# 
#     cd4_de_gwas = cd4_de.loc[[x.lower() in gwas for  x in cd4_de.gene]]
# 
#     plot_data_genes = [x.lower() for x in plot_data.index]
# 
#     ix = cd4_de_gwas.loc[[x.lower() in plot_data_genes for x in cd4_de_gwas.gene]]
#     all_res.append(ix)
# 
# all_res = pd.concat(all_res, axis=0)
# 
# 
# all_res[['ident', 'gene', 'log2fc']].round(2)
