import os

import numpy as np
from sklearn.manifold import TSNE

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import torch
from scvi.dataset import CortexDataset, RetinaDataset

from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer

import pandas as pd
import feather

#gene_dataset = CortexDataset(save_path="cortex")

counts_file = snakemake.input['exp']
meta_file = snakemake.input['meta']

latent_file = snakemake.output['latent']
genes_file = snakemake.output['model_genes']
model_file = snakemake.output['model']
batch_encoding_file = snakemake.output['batch_encoding']


out_dir = os.path.dirname(latent_file)

if not os.path.exists(out_dir):
    os.mkdir(out_dir)


if counts_file.endswith('feather'):
    counts = feather.read_dataframe(counts_file).set_index('index')
else:
    counts = pd.read_table(counts_file, index_col=0)

counts = counts.astype("int32").T

meta = pd.read_table(meta_file, index_col=0)

meta = meta.loc[counts.index]

try:
    has_batch = snakemake.params['has_batch']
except AttributeError:
    has_batch = True

try:
    batch_var = snakemake.params['batch_var']
except AttributeError:
    batch_var = 'Batch'

if has_batch:
    batch = meta[batch_var].astype('category').cat.codes.values
else:
    batch = np.zeros(meta.shape[0])

# Save batch encoding for later loading
pd.DataFrame({'Batch': batch}, index=meta.index) \
    .to_csv(batch_encoding_file, sep="\t")

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
    COMPONENTS = int(snakemake.params['components'])
except AttributeError:
    COMPONENTS = 10  # Latent component count

try:
    LAYERS = int(snakemake.params['layers'])
except AttributeError:
    LAYERS = 1  # number of hidden layers

try:
    zinb = bool(snakemake.params['zinb'])
    if zinb:
        reconstruction_loss = "zinb"
    else:
        reconstruction_loss = "nb"
except AttributeError:
        reconstruction_loss = "zinb"

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
    gene_passes = filters.filter_genes_fano(scaled_counts.T,
                                            num_mad=filterParam2)
    counts = counts.loc[:, gene_passes]

# Load data into scvi
from scvi.dataset.dataset import GeneExpressionDataset

batch = batch.reshape((-1, 1)).astype('int64')
cvals = counts.values.astype('int64')
zz = GeneExpressionDataset.get_attributes_from_matrix(cvals, batch_indices=batch)

# zz[0]: int64, ndarray, genes x cells
# zz[1]: float32, ndarray, cells x 1
# zz[2]: float32, ndarray, cells x 1
# zz[3]: int64, ndarray, cells x 1

dataset = GeneExpressionDataset(*zz, gene_names=counts.columns)

n_epochs = 400
lr = 1e-3
use_batches = True
use_cuda = torch.cuda.is_available()

torch.set_num_threads(20)
# However, need to set MKL_NUM_THREADS too

print("Running scVI with {} components, on {} batches and {} genes".format(
    COMPONENTS, dataset.n_batches, cvals.shape[1])
)

vae = VAE(dataset.nb_genes, n_batch=dataset.n_batches * use_batches,
          n_latent=COMPONENTS, n_layers=LAYERS,
          reconstruction_loss=reconstruction_loss)
trainer = UnsupervisedTrainer(vae,
                              dataset,
                              train_size=0.75,
                              use_cuda=use_cuda,
                              frequency=5)
trainer.train(n_epochs=n_epochs, lr=lr)


# Plot training result

plt.figure()
ll_train_set = trainer.history["ll_train_set"]
ll_test_set = trainer.history["ll_test_set"]
x = np.linspace(0, n_epochs, (len(ll_train_set)))
plt.plot(x, ll_train_set)
plt.plot(x, ll_test_set)
ymax = np.percentile(ll_test_set, 95)
ymin = np.min(ll_train_set) - .5*(ymax - np.min(ll_train_set))
plt.ylim(ymin, ymax)
plt.savefig(os.path.join(out_dir, 'training.png'))

tData = pd.DataFrame({
    'Training_Loss': ll_train_set,
    'Test_Loss': ll_test_set,
})

tData.to_csv(os.path.join(out_dir, 'training.txt'), sep="\t")


# Get latent space
latent = trainer.get_all_latent_and_imputed_values()["latent"]

# Save Results
latent = pd.DataFrame(latent, index=counts.index)
latent.to_csv(latent_file, sep="\t", compression="gzip")

torch.save(trainer.model, model_file)

genes = counts.columns.to_frame()
genes.to_csv(genes_file, sep="\t", index=False, header=False)
