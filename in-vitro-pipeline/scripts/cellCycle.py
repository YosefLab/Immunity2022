import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gene_enrich
import os
import subprocess

repo_dir = os.path.dirname(subprocess.check_output(['git', 'rev-parse', '--git-dir'])).decode()


# tpm_in ="data_filtered/tpm.txt.gz" 

tpm_in = snakemake.input['tpm']

cc_out = snakemake.output['cc']

out_dir = os.path.dirname(cc_out)
os.makedirs(out_dir, exist_ok=True)

tpm = pd.read_table(tpm_in, index_col=0)

gs = gene_enrich.load_gene_set_gmt(
    os.path.join(repo_dir, "Signatures/CellCycle/cell_cycle_Tirosh.gmt")
)
g1s = [x for x in gs if x.name == 'G1/S'][0]
g2m = [x for x in gs if x.name == 'G2/M'][0]

cc_genes = pd.Index(gs[0].genes | gs[1].genes) & tpm.index
cc_exp = tpm.loc[cc_genes]

cc_exp_norm = cc_exp.subtract(cc_exp.mean(axis=1), axis=0) \
                    .divide(cc_exp.std(axis=1), axis=0)


from sklearn.decomposition import PCA

model = PCA(n_components=20)
model.fit(cc_exp_norm.values.T)
pca = model.transform(cc_exp_norm.values.T)

pca = pd.DataFrame(pca, index=cc_exp_norm.columns)
pca = pca.iloc[:, 0:2]
pca.columns = ['CC1', 'CC2']
pca.to_csv(cc_out, sep='\t', compression='gzip')

plt.figure()
plt.bar(x=np.arange(20)+1, height=model.explained_variance_)
plt.ylabel('Component Scores')
plt.xlabel('Principal Components')
plt.savefig(os.path.join(out_dir, 'Components.png'))

row_colors = pd.DataFrame({
    'G1/S': ['#000000' if x in g1s.genes else '#FFFFFF' for x in cc_exp_norm.index],
    'G2/M': ['#000000' if x in g2m.genes else '#FFFFFF' for x in cc_exp_norm.index],
}, index=cc_exp_norm.index)

sns.clustermap(cc_exp_norm, vmin=-2, vmax=2, cmap="RdBu_r", row_colors=row_colors)
plt.savefig(os.path.join(out_dir, 'CC_Heatmap.png'))


# Plot cycling proportions
g1_exp = cc_exp_norm.loc[pd.Index(g1s.genes) & cc_exp_norm.index]
g2_exp = cc_exp_norm.loc[pd.Index(g2m.genes) & cc_exp_norm.index]

g1_score = g1_exp.mean(axis=0)
g2_score = g2_exp.mean(axis=0)


g1_thresh = -0.25
g2_thresh = 0.20

cycling = ((g1_score > g1_thresh) | (g2_score > g2_thresh)).rename('Cycling')
isg1 = ((g1_score > g1_thresh)).rename('G1')
isg2 = ((g2_score > g2_thresh)).rename('G2')
isg0 = (~cycling).rename('G0')


meta = pd.read_table("data_filtered/meta.txt.gz", index_col=0)
meta['StimGeno'] = [
    '{}:{}'.format(x, y) for x, y in zip(meta.Stimulus, meta.Genotype)
]
meta = meta.join(cycling).join(isg1).join(isg2).join(isg0)
cycle_count = meta.groupby(['StimGeno', 'GFP'])['Cycling'].sum()
prop = meta.groupby(['StimGeno', 'GFP']).size()
cc = (cycle_count/prop*100).reset_index()

plt.figure(figsize=(7, 4))
sns.barplot(x='StimGeno', hue='GFP', y=0, data=cc)
plt.ylabel('Cycling (%)')
plt.savefig(os.path.join(out_dir, 'CC_Proportions.png'))
