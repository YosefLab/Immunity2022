import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import json

import scanpy.api as sc

exp_file = snakemake.input['exp']
meta_file = snakemake.input['meta']
latent_file = snakemake.input['latent']
tsne_file = snakemake.input['tsne']

cluster_colors_file = snakemake.output['cluster_colors']
cluster_plot_file = snakemake.output['cluster_plot']
clusters_file = snakemake.output['clusters']

out_dir = os.path.dirname(clusters_file)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

data = pd.read_table(latent_file, index_col=0)


adata = sc.AnnData(data.values,
                   obs={'smp_names': data.index.tolist()},
                   var={'var_names': data.columns.tolist()})



try:
    n_neighbors = snakemake.params['n_neighbors']
except AttributeError:
    n_neighbors = 10

try:
    resolution = snakemake.params['resolution']
except AttributeError:
    resolution = 0.5

try:
    variable_name = snakemake.params['variable_name']
except AttributeError:
    variable_name = 'Cluster'

print('Running Louvain clustering with n_neighbors={} and resolution={}'.format(
    n_neighbors, resolution))

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)

sc.tl.louvain(adata, resolution=resolution)


clusters = adata.obs['louvain']
clusters = ['c' + x for x in clusters]
sample_names = adata.obs_names
clusters = pd.Series(clusters, index=sample_names, name=variable_name)
clusters.head()


# ------------------------------------------
# label some clusters based off expression
# ------------------------------------------

# edata = pd.read_table(exp_file, index_col=0)
# edata = edata.subtract(edata.mean(axis=1), axis=0)
# denom = edata.std(axis=1)
# denom[denom == 0] = 1.0
# edata = edata.divide(denom, axis=0)
# edata = edata.T.join(clusters)
#
#
# new_labels = {}
#
#
# def add_label(genes, new_label):
#
#     if isinstance(genes, str):
#         genes = [genes]
#
#     plt.figure()
#     vals = edata[genes].sum(axis=1)
#     old_label = pd.concat((vals, edata[variable_name]), axis=1).groupby(
#         variable_name).mean()[0].idxmax()
#     sns.boxplot(x=edata[variable_name], y=vals)
#     plt.xlabel(', '.join(genes))
#     plt.title(old_label + "-> " + new_label)
#     new_labels[old_label] = new_label
#
#
# add_label('RORC', 'Th17_like')
# add_label(['ITGB1', 'CCR2', 'S100A4', 'GATA3'], 'r0')
# add_label(['CCR7'], 'CCR7')
# add_label(['MKI67', 'BUB1'], 'r1')
# add_label(['CCL5', 'GZMA'], 'CCL5_GZMA')
# add_label(['CCL5', 'CCR5', 'GZMK'], 'CCL5_CCR5')
# add_label(['CCL4'], 'CCL4')


# give cluster's specific colors (to be used later)
if len(clusters.unique()) > 10:
    cmap = plt.get_cmap('tab20')
else:
    cmap = plt.get_cmap('tab10')

cluster_colors = {}
for i, cluster in enumerate(clusters.unique()):
    cluster_colors[cluster] = cmap.colors[i % len(cmap.colors)]


# Load the tSNE coordinates for plotting
plt.style.use('default')
tsne = pd.read_table(tsne_file, index_col=0)
tsne.columns = ['tsne1', 'tsne2']
tsne = tsne.join(clusters)

ms = 3
if tsne.shape[0] > 10000:
    ms = 2
if tsne.shape[0] > 30000:
    ms = 1

ff = plt.figure(figsize=(10, 5))
for clust in tsne[variable_name].unique():
    samples = tsne.loc[tsne[variable_name] == clust]
    color = cluster_colors[clust]
    plt.plot(samples['tsne1'], samples['tsne2'],
             'o', ms=ms, label=clust, color=color)

plt.subplots_adjust(right=.6)
plt.legend(loc='center left', bbox_to_anchor=(1, .5), markerscale=3/ms)
plt.savefig(cluster_plot_file)


# Save clusters to disk
clusters.to_frame().to_csv(clusters_file, sep='\t')

# Save cluster colors to disk
json.dump(cluster_colors, open(cluster_colors_file, "w"))
