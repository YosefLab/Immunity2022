# In[1]:

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import feather
from natsort import natsorted
from matplotlib.colors import LinearSegmentedColormap, Normalize
plt.rcParams['svg.fonttype'] = 'none'
import common


def despine():
    ax = plt.gca()
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# %% Load Data


meta = pd.read_table("analysis/data/meta.txt.gz", index_col=0)
clusters = pd.read_table("analysis/cluster_together/clusters.txt", index_col=0)

meta = meta.join(clusters)

# %% Get Proportions for each cluster

msub_ctrl = meta.loc[meta.Genotype == "ctrl"]
cluster_props_ctrl = msub_ctrl.Cluster.value_counts() / msub_ctrl.shape[0]

msub_ko = meta.loc[meta.Genotype == "KO"]
cluster_props_ko = msub_ko.Cluster.value_counts() / msub_ko.shape[0]

cluster_props = meta.Cluster.value_counts() / meta.shape[0] * 100

# relative change in proportion, ctrl vs ko
cluster_delta = (cluster_props_ctrl / cluster_props_ko)

# %% Plot results

cluster_order = cluster_delta.sort_values(ascending=False).index

cluster_props = cluster_props[cluster_order]
cluster_delta = cluster_delta[cluster_order]

fig, axs = plt.subplots(2, 1, figsize=(5, 5),
                        gridspec_kw=dict(height_ratios=[.5, 1])
                        )

plt.sca(axs[0])

axs[0].set_axisbelow(True)

plt.bar(
    x=[common.name_map_lpl_clusters[x] for x in cluster_props.index],
    height=cluster_props.values,
    color=[common.cmap_lpl_clusters[x] for x in cluster_props.index]
)
plt.xticks(rotation=45, size=7)

plt.text(-3, 10, 'Cluster\nProportion (%)',
         horizontalalignment='right', verticalalignment='center',
         multialignment='center')

despine()

plt.grid(axis='y', ls=(0, (5, 5)), lw=.5, color='#aaaaaa')
plt.yticks([4, 8, 12])

plt.sca(axs[1])

axs[1].set_axisbelow(True)

cmap = LinearSegmentedColormap.from_list(
    "test", [common.cmap_genotype['KO'], common.cmap_genotype['ctrl']]
)
norm = Normalize(vmin=-1, vmax=1, clip=True)
colors = [cmap(norm(np.log2(x))) for x in cluster_delta.values]

def tx_signed_height(x):
    if x >= 1:
        return x-1
    else:
        return 1-x

def tx_signed_bottom(x):
    if x >= 1:
        return 1
    else:
        return x


plt.bar(
    x=[common.name_map_lpl_clusters[x] for x in cluster_delta.index],
    height=[tx_signed_height(x) for x in cluster_delta.values],
    bottom=[tx_signed_bottom(x) for x in cluster_delta.values],
    color=colors
)
plt.yscale('log')
plt.yticks(
    ticks=[1/3, .5, 2/3, 1, 1.5, 2, 3],
    labels=['3x', '2x', '1.5x', '1:1', '1.5x', '2x', '3x']
)
plt.ylim(1/4, 4)

axs[1].tick_params(
    axis='y', which='minor', labelleft=False, left=False
)

axs[1].tick_params(
    axis='x', which='major', labelbottom=True, bottom=False
)

plt.text(-3, 2, 'Higher in\n$IL23R^{wt/eGFP}$',
         horizontalalignment='right', verticalalignment='center',
         multialignment='center')

plt.text(-3, 1/2, 'Higher in\n$IL23R^{eGFP/eGFP}$',
         horizontalalignment='right', verticalalignment='center',
         multialignment='center')


despine()
axs[1].spines['bottom'].set_visible(False)

# Add the x labels to the bars
plt.xticks(rotation=45, size=7)

# Alternate x ticks

# plt.xticks([])
# for i in range(len(cluster_delta)):
#     cl = common.name_map_lpl_clusters[cluster_delta.index[i]]
#     pos = cluster_delta.values[i] > 1
#     x = i-.3 if pos else i+.3
#     y = 1/1.05 if pos else 1.05
#     plt.text(x, y, cl, rotation=45, size=7,
#              verticalalignment='top' if pos else 'bottom',
#              horizontalalignment='center')
# 
# 
# plt.xlim(-1.5, len(cluster_delta))
plt.grid(axis='y', ls=(0, (5, 5)), lw=.5, color='#aaaaaa')
plt.grid(axis='x', ls='-', lw=.5, color='#dddddd')

axs[1]

plt.subplots_adjust(left=0.4, hspace=0.4)

#plt.show()
plt.savefig('ClusterProportions.svg')
