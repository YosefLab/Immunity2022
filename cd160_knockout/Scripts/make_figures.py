#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# %% Define some colors and styles

COOL_RED = sns.color_palette('coolwarm')[5]
COOL_BLUE = sns.color_palette('coolwarm')[0]

sns.set_style('whitegrid')
sns.set_context('notebook')
plt.rcParams['figure.dpi'] = 100
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['legend.framealpha'] = 1.0
plt.rcParams['patch.edgecolor'] = 'none'
plt.rcParams['svg.fonttype'] = 'none'

# %% Load the data

tpm = pd.read_table("rnaseq/tpm.txt", index_col=0)
meta = pd.read_table("rnaseq/meta.txt", index_col=0)

# %% Drop S1 - outlier sample

valid_samples = meta.drop('S1', axis=0).index
tpm = tpm.loc[:, valid_samples]
meta = meta.loc[valid_samples]


# %% Load DE Results

de_results = pd.read_table("de_genotype.txt", index_col=0)
de_results = de_results.drop('TDTOMATO', axis=0)


# %% Plot up-down genes

to_plot = de_results.loc[
    lambda x: (
        (x['FDR'] < .05) &
        (x['logFC'].abs() > .25)
    )
]


n_up = to_plot.loc[lambda x: x['logFC'] > 0].shape[0]
n_dn = to_plot.loc[lambda x: x['logFC'] < 0].shape[0]


fig = plt.figure(figsize=(4, 5))

plt.bar(
    ['Downregulated', 'Upregulated'],
    [n_dn, n_up],
    width=0.6,
    color=[COOL_BLUE, COOL_RED]
)

plt.grid(False, axis='x')
plt.gca().set_axisbelow(True)
plt.xlim(-.5, 1.5)

plt.ylabel('# of Genes')

plt.title('$CD160^{-/-}$ vs Ctrl\nDifferentially Expressed Genes\n(FDR < .05, |$log_2FC$| > .25)')

plt.subplots_adjust(bottom=0.3, left=0.2)
plt.savefig(
    'Figures/up_down_genes.svg', dpi=200
)
plt.close()

# %%

np.log(tpm.loc['SELL'] + 1)


# %% Volcano comparison with c9 1vsAll

cluster_1vall_results_file = "../10xAll/combined_LPL/cluster_together/cluster_de_1vAll_edgeR//cluster_de_results.txt.gz"
cluster_1vall_results = pd.read_table(cluster_1vall_results_file)

c9_up_genes = cluster_1vall_results.loc[
    lambda x: (
        (x['cluster'] == 'c9') &
        (x['logFC'] > .25) &
        (x['FDR'] < .05)
    )
].GeneSymbol.str.upper().tolist()
c9_dn_genes = cluster_1vall_results.loc[
    lambda x: (
        (x['cluster'] == 'c9') &
        (x['logFC'] < -.25) &
        (x['FDR'] < .05)
    )
].GeneSymbol.str.upper().tolist()
sig_genes_cd160 = de_results.loc[
    lambda x: (
        (x['logFC'].abs() > .25) &
        (x['FDR'] < .05)
    )
].index.str.upper().tolist()


de_results['Group'] = 'n.s.'
de_results.loc[
    lambda x: (
        x.index.isin(sig_genes_cd160)
    ), 'Group'
] = 'DE (CD160)'
de_results.loc[
    lambda x: (
        x.index.isin(c9_up_genes) & x.index.isin(sig_genes_cd160)
    ), 'Group'
] = 'DE (CD160) and upregulated in C9'
de_results.loc[
    lambda x: (
        x.index.isin(c9_dn_genes) & x.index.isin(sig_genes_cd160)
    ), 'Group'
] = 'DE (CD160) and downregulated in C9'


cmap = {
    'n.s.': '#bbbbbb55',
    'DE (CD160)': '#bbbbbb55',
    'DE (CD160) and upregulated in C9': COOL_RED,
    'DE (CD160) and downregulated in C9': COOL_BLUE,
}
sizemap = {
    'n.s.': 9,
    'DE (CD160)': 9,
    'DE (CD160) and upregulated in C9': 15,
    'DE (CD160) and downregulated in C9': 15,
}

sns.color_palette('coolwarm')[5]
de_results['Color'] = de_results['Group'].map(cmap)
de_results['Size'] = de_results['Group'].map(sizemap)


# %% Plot Volcano


plt.figure(figsize=(6, 6))
plt.scatter(
    x=de_results['logFC'],
    y=np.log10(de_results['FDR'])*-1,
    s=de_results['Size'],
    edgecolors='none',
    c=de_results['Color'],
    rasterized=True
)
plt.ylim(0, plt.ylim()[1])
plt.ylabel('Adjusted P-Value ($-log_{10}$)')
plt.xlabel('Fold-Change ($log_2$)\nCD160 KO vs. Ctrl')
plt.grid(lw=0.5)

from matplotlib.lines import Line2D
labels = ['DE (CD160) and upregulated in C9', 'DE (CD160) and downregulated in C9']
handles = [
    Line2D([0], [0], color=cmap[x], markersize=sizemap[x]**.5, marker='o', lw=0)
    for x in labels
]
plt.legend(handles, labels, loc='lower left', bbox_to_anchor=[0, 1], frameon=False)
plt.savefig('Figures/volcano_c9_sig.svg', dpi=200)
plt.close()


# %% Enrichments

gsea_res = pd.read_table("gsea/res_table.txt", index_col='Term')

# **Notable Enrichments**
enrichments_to_plot = [
    'HALLMARK_INFLAMMATORY_RESPONSE',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE',
    'TYROBP_Causal_Network',
    'toll-like receptor signaling pathway',
    'cell adhesion',
    'chemotaxis',
    'positive regulation of cell migration',
    'positive regulation of Wnt signaling pathway',
    'cytokine-mediated signaling pathway'
]

y_axis = [x.title().replace('_', ' ') for x in enrichments_to_plot]
y_axis = [x.replace('Tyrobp', 'TYROBP') for x in y_axis]
x_axis = np.log10(gsea_res.loc[enrichments_to_plot]['fdr'])*-1

fig = plt.figure(figsize=(7, 7))

plt.barh(
    y_axis, x_axis, height=0.5, color='#1ea1c2'
)

plt.subplots_adjust(left=0.5)
plt.yticks(size=8)
plt.xlabel('Adjusted P-Value ($-log_{10}$)', size=10)
plt.grid(False, axis='y')

plt.savefig('Figures/enrichments.svg', dpi=200)
plt.close()


# %% Signature scores



def get_sigs(df):
    sigs = {}
    for cluster in df['cluster'].unique():
        sub = df.loc[
            lambda x: (
                (x['cluster'] == cluster) &
                (x['FDR'] < .05) &
                (x['logFC'].abs() > .25)
            )
        ]

        sig_up_genes = sub.loc[lambda x: x['logFC'] > 0].GeneSymbol.tolist()
        sig_dn_genes = sub.loc[lambda x: x['logFC'] < 0].GeneSymbol.tolist()

        sig_up_genes = [x.upper() for x in sig_up_genes]
        sig_dn_genes = [x.upper() for x in sig_dn_genes]

        sig = {x: 1 for x in sig_up_genes}
        sig.update({x: -1 for x in sig_dn_genes})
        sig = pd.Series(sig)
        sigs[cluster] = sig

    return sigs

sigs_1vall = get_sigs(cluster_1vall_results)
sigs_1vall = {k+'-1vall': v for k, v in sigs_1vall.items()}

# In[121]:

log_tpm = np.log(tpm + 1)
log_tpm_norm = log_tpm.subtract(
    log_tpm.mean(axis=0), axis=1).divide(
    log_tpm.std(axis=0), axis=1
)
log_tpm_norm = log_tpm_norm.drop('CD160', axis=0)  # Can't use this to test as we knocked it out specifically


# In[122]:


def score_sig(sig_name, sig_values):
    sig_up_genes = sig_values.index[sig_values > 0] & log_tpm_norm.index
    sig_dn_genes = sig_values.index[sig_values < 0] & log_tpm_norm.index

    sig_up_score = log_tpm_norm.loc[sig_up_genes].mean(axis=0)
    sig_dn_score = log_tpm_norm.loc[sig_dn_genes].mean(axis=0)
    sig_score = sig_up_score - sig_dn_score
    sig_score.name = sig_name
    return sig_score

sig_scores = [score_sig(name, vals) for name, vals in sigs_1vall.items()]
sig_scores = pd.concat(sig_scores, axis=1)

plot_data = sig_scores.join(meta)

from scipy.stats import ttest_ind


pvals = []
sigs_to_test = ['c2-1vall', 'c7-1vall', 'c9-1vall']
for sig in sigs_to_test:
    a = plot_data.loc[
        lambda x: x['Genotype'] == 'ctrl', sig
    ]
    b = plot_data.loc[
        lambda x: x['Genotype'] == 'KO', sig
    ]

    ttest_result = ttest_ind(a.values, b.values)
    pvals.append(ttest_result.pvalue)

from statsmodels.stats.multitest import multipletests
pvals_corrected = multipletests(pvals, method='fdr_bh')[1]
pvals_corrected = pd.Series(pvals_corrected, index=sigs_to_test)


# %%

colors = plt.get_cmap('Dark2').colors
genotype_map = {
    'ctrl': 'Ctrl',
    'KO': '$CD160^{-/-}$',
}
plot_data['GenotypeNice'] = plot_data['Genotype'].map(genotype_map)
plt.figure(figsize=(5, 5))
sns.boxplot(
    data=plot_data,
    x='GenotypeNice',
    y='c9-1vall',
    color=colors[0]
)
plt.xlabel('')
plt.ylabel('C9 vs. All\nSignature Score')
plt.subplots_adjust(left=0.3, bottom=0.2, right=0.85)
plt.savefig(
    'Figures/sig_scores_1vall_c9.svg', dpi=200
)
plt.close()

# %% Heatmap
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list

genes_to_plot = de_results.loc[lambda x: (
    (x['logFC'].abs() > .5) &
    (x['FDR'] < .1)
)].index


logtpm = np.log(tpm + .1)
lognormtpm = logtpm.subtract(logtpm.mean(axis=1), axis=0)

plot_data = lognormtpm.loc[genes_to_plot]

# Determine gene order
pos_genes = de_results.loc[genes_to_plot].loc[lambda x: x['logFC'] > 0].index
neg_genes = de_results.loc[genes_to_plot].loc[lambda x: x['logFC'] < 0].index


order_pos = pos_genes[leaves_list(linkage(plot_data.loc[pos_genes]))]
order_neg = neg_genes[leaves_list(linkage(plot_data.loc[neg_genes]))]
order = order_pos.tolist() + order_neg.tolist()

plot_data = plot_data.loc[order]

# %% plot it

plot_data.index = [x.capitalize() for x in plot_data.index]

genes_to_highlight = {
    "Saa3", "Steap4", "Cd1d2", "Acod1",
    "Sell", "Ccl1", "Il10ra", "Ets1",
    "Cd160", "Ddc",
}

# cmap = LinearSegmentedColormap.from_list("purple", ["#FF00FF", "#000000", "#FFFF00"])
cmap = LinearSegmentedColormap.from_list("coolbluewhitered", [COOL_BLUE, "#FFFFFF", COOL_RED])
cm = sns.clustermap(
    plot_data,
    vmin=-1.5, vmax=1.5, cmap=cmap,
    col_cluster=False, row_cluster=False,
    yticklabels=True,
    xticklabels=False,
    rasterized=True
)

plt.sca(cm.ax_heatmap)
locs, labels = plt.yticks()

new_locs, new_labels = [], []

for loc, label in zip(locs, labels):
    if label.get_text() in genes_to_highlight:
        new_locs.append(loc)
        new_labels.append(label)

plt.yticks(new_locs, new_labels, size=9)

plt.savefig("Figures/gene_expression_heatmap.svg", dpi=200)
