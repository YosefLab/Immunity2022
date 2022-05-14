import os
import pandas as pd
import scipy.io
import scipy.sparse
import h5py
import numpy as np

ens_file = snakemake.output['ens']
sym_file = snakemake.output['sym']
filtered_file = snakemake.output['filtered']
scaled_file = snakemake.output['scaled']
genes_file = snakemake.output['genes']

# make the output directory
os.makedirs(os.path.dirname(ens_file), exist_ok=True)

cellranger_output = snakemake.input['cellranger_out']

f = h5py.File(cellranger_output)

barcodes = f['mm10/barcodes'][()].astype(np.str_)
barcodes = pd.Series(barcodes, name='barcode')

gs = pd.Series(f['mm10/gene_names'][()].astype(np.str_), name='GeneSymbol')
gi = pd.Series(f['mm10/genes'][()].astype(np.str_), name='geneID')
genes = pd.concat((gi, gs), axis=1).set_index('geneID')
genes['GeneSymbol'] = genes['GeneSymbol'].str.upper()

genes.to_csv(genes_file, sep="\t")

data = f['mm10/data'][()]
indices = f['mm10/indices'][()]
indptr = f['mm10/indptr'][()]
shape = f['mm10/shape'][()]
gx = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
gx = gx.toarray()

gx = pd.DataFrame(gx, index=genes.index, columns=barcodes)

gx.to_csv(ens_file, sep='\t', compression='gzip')

gx_sym = gx.join(genes).groupby('GeneSymbol').sum()

gx_sym.to_csv(sym_file, sep='\t', compression='gzip')

# filter

# For now, just remove genes with less than 10 counts in total
valid_genes = gx_sym.sum(axis=1) > 10

gx_filtered = gx_sym.loc[valid_genes]

gx_filtered.to_csv(filtered_file, sep='\t', compression='gzip')


# scale

counts_per_cell = gx_filtered.sum(axis=0)

gx_norm = gx_filtered.divide(counts_per_cell, axis=1)*4000

gx_norm.to_csv(scaled_file, sep='\t', compression='gzip')
