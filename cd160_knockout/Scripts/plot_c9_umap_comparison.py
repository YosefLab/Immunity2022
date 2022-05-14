import numpy as np
import pandas as pd
import feather
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('whitegrid')
sns.set_context('notebook')
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['legend.framealpha'] = 1.0
plt.rcParams['patch.edgecolor'] = 'none'
plt.rcParams['svg.fonttype'] = 'none'

# %% Load the CD160 results and define the signature

cd160_de = pd.read_table("de_genotype.txt", index_col=0)
cd160_de.index = [x.capitalize() for x in cd160_de.index]

up_genes = cd160_de.loc[lambda x: (
    (x['FDR'] < .1) &
    (x['logFC'] > .5)
)].index

dn_genes = cd160_de.loc[lambda x: (
    (x['FDR'] < .1) &
    (x['logFC'] < .5)
)].index

dn_genes = dn_genes.difference(['Cd160'])
print(len(up_genes), len(dn_genes))


# %% Load the gene expression data

scaled_exp = feather.read_dataframe(
    "../10xAll/combined_LPL/data/expression_scaled.feather"
)
scaled_exp.index = scaled_exp['index']
scaled_exp.index = [x.capitalize() for x in scaled_exp.index]
scaled_exp = scaled_exp.drop('index', axis=1)

clusters = pd.read_table(
    "../10xAll/combined_LPL/cluster_together/clusters.txt", index_col=0
)
clusters = clusters.loc[scaled_exp.columns]
cluster_colors = json.load(open("../10xAll/combined_LPL/cluster_together/cluster_colors.json"))

umap = pd.read_table(
    "../10xAll/combined_LPL/umap/umap.txt", index_col=0
)
umap = umap.loc[scaled_exp.columns]

# %% Evaluate a signature score

logexp = np.log(scaled_exp+1)

up_score = logexp.loc[[x in up_genes for x in logexp.index], :].mean(axis=0)
dn_score = logexp.loc[[x in dn_genes for x in logexp.index], :].mean(axis=0)

sig_score = (up_score - dn_score)

# %% What are scores vs cluster?

scores_c9 = sig_score.loc[clusters.Cluster == 'c9']
scores_other = sig_score.loc[clusters.Cluster != 'c9']

scores_df = clusters.copy()
scores_df['sig_score'] = sig_score
scores_df['cluster_group'] = 'Other'
scores_df.loc[lambda x: x['Cluster'] == 'c9', 'cluster_group'] = 'c9'
scores_df = scores_df.join(umap)

from scipy.stats import ttest_ind
ttest_ind(scores_c9, scores_other)

# %%
from matplotlib.gridspec import GridSpec

from matplotlib.lines import Line2D

scores_df['color'] = [cluster_colors['c9'] if c == 'c9' else '#CCCCCC' for c in scores_df['Cluster']]
scores_df['size'] = [2 if c == 'c9' else .5 for c in scores_df['Cluster']]


fig = plt.figure(figsize=(9, 9))
ax_cluster = fig.add_axes([.05, .6, .3, .3])
ax_sig = fig.add_axes([.05, .1, .5, .5])
ax_box = fig.add_axes([.7, .2, .28, .28])

plt.sca(ax_cluster)

plt.scatter(
    x=scores_df['umap1'], y=scores_df['umap2'],
    c=scores_df['color'], s=scores_df['size'],
    rasterized=True
)
plt.xticks([])
plt.yticks([])
sns.despine(left=True, bottom=True, ax=plt.gca())

handles = [
    Line2D([0], [0], marker='o', markerfacecolor=cluster_colors['c9'], markersize=5, linewidth=0, markeredgewidth=0),
    Line2D([0], [0], marker='o', markerfacecolor='#CCCCCC', markersize=5, linewidth=0, markeredgewidth=0)
]
labels = ['C9', 'Other']

fig.legend(handles=handles, labels=labels, loc='lower left', bbox_to_anchor=[.35, .65], frameon=False, title='Cluster')

plt.sca(ax_sig)

plt.scatter(
    x=scores_df['umap1'], y=scores_df['umap2'], c=scores_df['sig_score'], s=1, vmin=-.18, vmax=-.03, cmap='viridis',
    rasterized=True
)
sns.despine(left=True, bottom=True, ax=plt.gca())
plt.xticks([])
plt.yticks([])
cax = fig.add_axes([.07, .1, .02, .1])
plt.colorbar(cax=cax)
cax.yaxis.set_label_position("left")
cax.set_ylabel('$CD160^{-/-}$ vs WT\nSignature Score')

plt.sca(ax_box)
sns.boxplot(data=scores_df, x='cluster_group', y='sig_score', order=['c9', 'Other'], palette=[cluster_colors['c9'], '#CCCCCC'])
plt.xlabel('Cluster')
plt.ylabel('$CD160^{-/-}$ vs WT\nSignature Score')
sns.despine(ax=ax_box, left=True, right=False)
ax_box.yaxis.set_label_position("right")
ax_box.yaxis.tick_right()
plt.savefig('Figures/lpl_umap.svg', dpi=200)
