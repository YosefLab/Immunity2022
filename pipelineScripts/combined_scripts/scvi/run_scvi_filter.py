# Same as run_scvi only more parameters and options for 
# filtering the genes
import numpy as np
import pandas as pd
import scvi
import feather
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


counts_file = snakemake.input['exp']
meta_file = snakemake.input['meta']
latent_file = snakemake.output['latent']
genes_file = snakemake.output['model_genes']


out_dir = os.path.dirname(latent_file)

if not os.path.exists(out_dir):
    os.mkdir(out_dir)


if counts_file.endswith('feather'):
    counts = feather.read_dataframe(counts_file).set_index('index')
else:
    counts = pd.read_table(counts_file, index_col=0)

counts = counts.astype("int32").T

meta = pd.read_table(meta_file, index_col=0)

try:
    has_batch = snakemake.params['has_batch']
except AttributeError:
    has_batch = True

if has_batch:
    batch = meta.Batch.astype('category').cat.codes.values
else:
    batch = np.zeros(meta.shape[0])

try:
    filterType = snakemake.params['filterType']
except AttributeError:
    filterType = 'Threshold'

try:
    filterParam1 = snakemake.params['filterParam1']
except AttributeError:
    filterParam1 = 100  # N_CELLS

try:
    filterParam2 = snakemake.params['filterParam2']
except AttributeError:
    filterParam2 = 2  # MAD for fano

try:
    apply_mmd = snakemake.params['apply_mmd']
except AttributeError:
    apply_mmd = False  # Whether to apply mmd

# do some filtering
assert filterType in ['Threshold', 'Fano']

if filterType == 'Threshold' or filterType == 'Fano':
    N_COUNTS = 1
    N_CELLS = filterParam1

    gene_passes = (counts >= N_COUNTS).sum(axis=0) >= N_CELLS

    counts = counts.loc[:, gene_passes]

if filterType == 'Fano':
    from bio_utils import filters
    umis = counts.sum(axis=1)
    scaled_counts = counts.divide(umis, axis=0)*umis.median()
    gene_passes = filters.filter_genes_fano(scaled_counts.T, num_mad=filterParam2)
    counts = counts.loc[:, gene_passes]

holdout = int(counts.shape[0]/10)
model = scvi.scVIModel(counts, batch, test_holdout=holdout, apply_mmd=apply_mmd)

train_history = model.fit(num_epochs=250)

# Plot the history
plt.figure()
plt.plot(train_history['epoch'], train_history['t_loss'],
         label='Training Loss')
plt.plot(train_history['epoch'], train_history['v_loss'],
         label='Validation Loss')
plt.legend(loc='best')
plt.title('scVI Training Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.savefig(os.path.join(out_dir, 'training.png'))

latent = model.generate_latent_space()

latent.to_csv(latent_file, sep="\t", compression="gzip")

model.save(out_dir)

genes = counts.columns.to_frame()
genes.to_csv(genes_file, sep="\t", index=False, header=False)
