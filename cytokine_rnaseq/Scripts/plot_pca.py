import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
sns.set_palette('Dark2')

meta = pd.read_table("rnaseq/meta.txt", index_col=0)

tpm = pd.read_table("rnaseq/tpm.txt", index_col=0)

# Remove outliers
outliers = {
    'S1', 'S12', 'S19', 'S20', 'S21', 'S33', 'S34'
}

valid_samples = meta.index[~meta.index.isin(outliers)]

tpm = tpm.loc[:, valid_samples]
meta = meta.loc[valid_samples]

# standardize
tpm_ctrl = tpm.loc[:, meta.Genotype == 'ctrl']
logtpm = np.log(tpm_ctrl + 1)
tpm_norm = logtpm
valid_genes = tpm.max(axis=1) > 10
tpm_norm = tpm_norm.loc[valid_genes]

model = PCA(n_components=10)
pca = model.fit_transform(tpm_norm.values.T)
pca = pd.DataFrame(
    pca, index=tpm_norm.columns,
    columns=['PC{}'.format(i+1) for i in range(pca.shape[1])]
)
components = pd.DataFrame(
    model.components_, index=pca.columns,
    columns=tpm_norm.index).T

plot_data = pca.join(meta)

fig = plt.figure(figsize=(7.5, 5.5))

condition_colors = plt.get_cmap('Dark2').colors
condition_colors = [condition_colors[i] for i in [0, 1, 2, 7]]

sns.scatterplot(
    x='PC1', y='PC2', hue='Condition', style='IL-23R',
    s=120, data=plot_data, palette=condition_colors, markers=['o', '^']
)
plt.legend(loc='upper left', bbox_to_anchor=[1, 1], frameon=False)

plt.subplots_adjust(wspace=0.3, hspace=0.3, right=0.65)
plt.xlabel('PC-1', size=12, labelpad=10)
plt.ylabel('PC-2', size=12)

plt.savefig('Figures/cytokine_pca.svg')

pca.to_csv("pca.csv")
