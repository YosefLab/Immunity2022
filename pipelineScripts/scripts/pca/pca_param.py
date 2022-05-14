import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import feather
from bio_utils import filters


MAX_COMPONENTS = 80
BG_FACTOR = 2

import sys
expression_file = snakemake.input['exp']
output_file = snakemake.output['pca']
output_file_image = snakemake.output['image']
num_pcs = snakemake.params['num_pcs']

out_dir = os.path.dirname(output_file)
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

if expression_file.endswith('feather'):
    expression = feather.read_dataframe(expression_file)
    expression = expression.set_index('index')
else:
    expression = pd.read_table(expression_file, index_col=0)

# Filter expression
gene_passes = filters.filter_genes_fano(expression)
expression = expression.loc[gene_passes]

expression = np.log2(expression + 1)


model = PCA(n_components=MAX_COMPONENTS, random_state=0)
model.fit(expression.values.T)


# Now, look at null data
null = expression.values.copy()
for i in range(null.shape[0]):
    np.random.shuffle(null[i, :])

null_model = PCA(n_components=MAX_COMPONENTS, random_state=0)
null_model.fit(null.T)


# First value less than = num PCs
# E.g., if index 21 is first value less than null, then 0-20
# are greater than null, and num_pcs = 21
# num_pcs = np.nonzero(
#     (model.explained_variance_ < null_model.explained_variance_*BG_FACTOR)
# )[0][0]

# This wasn't working well for this data.  Not sure why.
# Just going to hard-code 20 for now
# num_pcs = 20

# Diagnostic plot
xmin = 0
xmax = min(num_pcs*2, 200)
ymin = 0
ymax = max(model.explained_variance_.max(),
           null_model.explained_variance_.max())*1.2

plt.figure()
plt.plot(np.arange(len(model.explained_variance_))+1,
         model.explained_variance_, 'o-', label='Expression',
         ms=2)
plt.plot(np.arange(len(null_model.explained_variance_))+1,
         null_model.explained_variance_, 'o-', label='Null',
         ms=2)
# plt.plot(np.arange(len(null_model.explained_variance_))+1,
#          null_model.explained_variance_*BG_FACTOR, 'o-', label='+Null',
#          ms=2)

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

plt.vlines(num_pcs, ymin, ymax, linestyle='dashed', color='red')
plt.legend(loc='best')
plt.title('Significant PCs: {}\n Num Genes: {}'.format(num_pcs, expression.shape[0]))
plt.savefig(output_file_image)

pcs = model.fit_transform(expression.values.T)
pcs = pcs[:, 0:num_pcs]

pcs = pd.DataFrame(pcs, index=expression.columns)

pcs.to_csv(output_file, sep='\t')
