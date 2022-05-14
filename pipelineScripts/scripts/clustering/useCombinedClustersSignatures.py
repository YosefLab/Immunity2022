from gene_enrich import load_gene_set_gmt
import numpy as np
import pandas as pd
import shutil


# expression_file = "expression/expression_scaled.txt.gz"
# signature_file = "../combined_all/cluster_together/cluster_signatures.gmt"
# colors_file = "../combined_all/cluster_together/cluster_colors.json"
# 
# cluster_out = "clustering/clusters.txt"
# cluster_colors_out = "clustering/cluster_colors.json"

expression_file = snakemake.input['exp']
signature_file = snakemake.input['signatures']
colors_file = snakemake.input['colors']

cluster_out = snakemake.output['clusters']
cluster_colors_out = snakemake.output['colors']

expression = pd.read_table(expression_file, index_col=0)
expression = np.log2(expression + 1)

sigs = load_gene_set_gmt(signature_file)

# Build the sig Matrix
sigMatrix = pd.DataFrame(0.0, index=expression.index,
                         columns=[x.name for x in sigs])

for sig in sigs:

    values = sig.values
    values = values.loc[values.index & expression.index]
    sigMatrix.loc[values.index, sig.name] = values

classify = sigMatrix.T.dot(expression)

proba = 1 / (1 + np.exp(classify*-1))

clusters = proba.idxmax(axis=0)
clusters = clusters.to_frame()
clusters.columns = ['Cluster']

clusters.to_csv(cluster_out, sep="\t")
shutil.copy(colors_file, cluster_colors_out)
