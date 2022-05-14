import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

plt.rcParams['svg.fonttype'] = 'none'
sns.set_palette('Dark2')

meta = pd.read_table("rnaseq/meta.txt", index_col=0)

tpm = pd.read_table("rnaseq/tpm.txt", index_col=0)

meta['FullName'] = [
    '{}-{}-{}'.format(x, y, z) for x, y, z in
    zip(meta['Condition'], meta['IL-23R'], meta['Genotype'])
]

# Remove outliers (low align % or read count)
outliers = {
    'S1', 'S12', 'S19', 'S20', 'S21', 'S33', 'S34'
}

valid_samples = meta.index[~meta.index.isin(outliers)]

tpm = tpm.loc[:, valid_samples]
meta = meta.loc[valid_samples]

# Load the ANOVA results
de_results = pd.read_table("de_anova.txt", index_col=0)
de_results = de_results.sort_values('PValue')
top_genes = de_results.sort_values('PValue').index[:200]

# Transform for plot
valid_genes = top_genes
valid_samples = meta.index[meta['Genotype'] == 'ctrl']
top_genes_log_tpm = np.log(tpm.loc[top_genes, valid_samples] + 1)
top_genes_log_tpm = top_genes_log_tpm.subtract(
    top_genes_log_tpm.mean(axis=1), axis=0
).divide(top_genes_log_tpm.std(axis=1), axis=0)


genes_to_plot = {
    "CCL20", "IL17F", "IFNG", "IL17A", "IL1R2", "NKG7",
    "CCL4", "GZMB", "IL9", "CCL3", "PENK", "PLAC8", "TGFBI",
    "IL12RB2", "IRF8", "LTA", "IGFBP4", "CD40LG", "CD83", "RORC",
    "CCR6", "BATF", "ITGA3", "EGFP_IRES_IL23R", "RORA", "IL10RA",
    "TBX21", "IL23R", "NFKBIA", "IL22", "SOCS2", "IRF7", "IFIT1",
    "STAT5A", "ID3", "CD74", "TCF7", "IL1R1", "IRF4", "IFNGR1",
    "TNFSF8", "IL2", "SELL", "IL18R1", "TNFRSF8", "CD28", "CCL5",
    "CTLA4", "CCR5", "LAG3", "IFITM1", "CCR7", "GPR18", "CD81",
    "FAS", "TNFRSF4", "CD86", "KLRC1", "GATA3", "CXCR3", "ITGB2",
    "IL12RB1", "TNF", "IL21", "IL21R", "IL7R", "IL18RAP", "CD44",
}

gfp_colors = {
    'positive': '#1B763B',
    'negative': '#BDBEC0'
}

cmap = LinearSegmentedColormap.from_list(
    "purple", ["#FF00FF", "#000000", "#FFFF00"]
)

col_order = meta.loc[top_genes_log_tpm.columns].sort_values(
    ['Condition', 'IL-23R'], ascending=[False, True]).index

top_genes_log_tpm = top_genes_log_tpm.loc[:, col_order]

col_colors = pd.DataFrame(
    {
        'GFP': meta.loc[top_genes_log_tpm.columns]['IL-23R'].map(gfp_colors)
    }, index=top_genes_log_tpm.columns
)

cm = sns.clustermap(
    top_genes_log_tpm,
    cmap=cmap, vmin=-2, vmax=2,
    cbar_kws={
        'ticks': [-2, 0, 2],
        'label': 'Standardized\n$Log_2$-Expression',
    },
    xticklabels=meta.loc[top_genes_log_tpm.columns].FullName,
    yticklabels=True,
    col_cluster=False,
    col_colors=col_colors,
    rasterized=True,
    figsize=(8, 10),
)

plt.sca(cm.ax_heatmap)
for text in cm.ax_heatmap.get_yticklabels():
    if text.get_text() not in genes_to_plot:
        text.set_visible(False)
    else:
        text.set_size(8)
        text.set_position((1, 200))

cm.ax_heatmap.tick_params(which='major', length=0)

# Create legend
labels = ['GFP+', 'GFP-']
handles = [
    Patch(facecolor=gfp_colors['positive']),
    Patch(facecolor=gfp_colors['negative']),
]
fig = plt.gcf()
fig.legend(
    handles, labels, bbox_to_anchor=[.3, .88], loc='upper left', frameon=False
)

plt.subplots_adjust(bottom=0.3)
plt.savefig('Figures/cytokine_heatmap.svg', dpi=300)
