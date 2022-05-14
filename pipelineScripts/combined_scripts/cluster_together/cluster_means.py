import os
import numpy as np
import pandas as pd
import feather


exp_file = snakemake.input['exp']
meta_file = snakemake.input['meta']
clusters_file = snakemake.input['clusters']

means_file = snakemake.output['out']


# Load data

expression = feather.read_dataframe(exp_file).set_index("index")
clusters = pd.read_table(clusters_file, index_col=0)

meta = pd.read_table(meta_file, index_col=0)

cluster_means = expression.T.join(clusters).groupby("Cluster").mean()
cluster_means['Genotype'] = 'All'
cluster_means = cluster_means.reset_index()

cg_means = expression.T.join(meta[['Genotype']].join(clusters)) \
    .groupby(["Genotype", "Cluster"]) \
    .mean()

cg_means = cg_means.reset_index()


cluster_means = cluster_means[cg_means.columns]

all_means = pd.concat(
    (cluster_means, cg_means),
    axis=0
)

all_means.to_csv(means_file, sep="\t", index=False)
