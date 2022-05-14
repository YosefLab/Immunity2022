import os
import pandas as pd
import scipy.io

ens_file = snakemake.output["ens"]
sym_file = snakemake.output["sym"]
filt_file = snakemake.output["filt"]
scale_file = snakemake.output["scale"]

os.makedirs(os.path.dirname(ens_file), exist_ok=True)

m_file = snakemake.input["matrix"]
b_file = snakemake.input["barcodes"]
g_file = snakemake.input["genes"]

gx = scipy.io.mmread(m_file)
gx = gx.toarray()

barcodes = pd.read_table(b_file, header=None)
barcodes = barcodes[0].tolist()

genes = pd.read_table(g_file, header=None)
gene_ids = genes[0].tolist()

gx = pd.DataFrame(gx, index=gene_ids, columns=barcodes)

gx.to_csv(ens_file, sep='\t', compression='gzip')


genes = genes.set_index(0)
genes.columns = ['GeneSymbol']
genes['GeneSymbol'] = genes['GeneSymbol'].str.upper()

gx_sym = gx.join(genes).groupby('GeneSymbol').sum()


gx_sym.to_csv(sym_file, sep='\t', compression='gzip')

# Filter
valid_genes = gx_sym.sum(axis=1) > 10

gx_filtered = gx_sym.loc[valid_genes]

gx_filtered.to_csv(filt_file, sep='\t', compression='gzip')

# Scale
counts_per_cell = gx_sym.sum(axis=0)

gx_norm = gx_filtered.divide(counts_per_cell, axis=1)*4000

gx_norm.to_csv(scale_file, sep='\t', compression='gzip')
