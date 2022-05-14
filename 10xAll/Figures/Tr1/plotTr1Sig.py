import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import common
import bio_utils
import seaborn as sns
import feather

plt.rcParams['svg.fonttype'] = 'none'

# %% Load Data

clusters = pd.read_table("analysis_lpl/cluster_together/clusters.txt", index_col=0)

sigScores = feather.read_dataframe("analysis_lpl/vision/sigScores.feather").set_index("index")


umap = pd.read_table("analysis_lpl/umap/umap.txt", index_col=0)

plot_data = clusters.join(sigScores[['Tr1_vs_other_CD4']]).join(umap)

plot_data['In7'] = 'Remainder'
plot_data.loc[plot_data.Cluster == 'c7', 'In7'] = 'LPL-7'

# %% Plot

fig, axs = plt.subplots(
    1, 2,
    figsize=(9, 5),
    gridspec_kw=dict(
        width_ratios=[1.8, 1],
        wspace=0.4
    )
)

cbar_ax = fig.add_axes(
    [.1, .2, 0.01, 0.15]
)

plt.sca(axs[0])

sc = plt.scatter(
    x=plot_data.umap1,
    y=plot_data.umap2,
    c=plot_data.Tr1_vs_other_CD4,
    s=3,
    cmap='viridis',
    vmin=-.2,
    vmax=.4,
    rasterized=True,
    edgecolors='none',
)


plt.xticks([])
plt.yticks([])
for sp in plt.gca().spines.values():
    sp.set_visible(False)

plt.colorbar(sc, cax=cbar_ax, ticks=[-.2, .4])
cbar_ax.set_yticklabels(['-0.2', ' 0.4'], size=8)
cbar_ax.set_ylabel('Tr1 Signature\nScore', size=9)
cbar_ax.yaxis.set_label_position('left')


plt.sca(axs[1])
colors = [
    common.cmap_lpl_clusters['c7'],
    '#CCCCCC'
]

sns.violinplot(data=plot_data, x='In7', y='Tr1_vs_other_CD4',
               order=['LPL-7', 'Remainder'], palette=colors)
plt.ylabel('Tr1 Signature Score')
plt.xlabel('')
plt.gca().set_axisbelow(True)
plt.grid(axis='y', color='#CCCCCC', linestyle=(0, (5, 5)))

# plt.show()
plt.savefig('tr1_signature.svg', dpi=300)
