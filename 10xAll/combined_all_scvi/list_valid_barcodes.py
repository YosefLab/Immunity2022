import pandas as pd

in_file = snakemake.input['clusters']
out_file = snakemake.output[0]

valid_clusters = {"c0", "c1", "c2", "c3", "c4", "c5", "c7"}


clusters = pd.read_table(in_file, index_col=0).Cluster

clusters_sub = clusters.loc[[x in valid_clusters for x in clusters]]

valid_barcodes = clusters_sub.index

with open(out_file, 'w') as fout:
    for bc in valid_barcodes:
        fout.write(bc)
        fout.write("\n")
