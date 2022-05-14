import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import common
from matplotlib.patches import Patch

plt.rcParams['svg.fonttype'] = 'none'

umap = pd.read_table("analysis/umap/umap.txt", index_col=0)
meta = pd.read_table("analysis/data/meta.txt.gz", index_col=0)


# %% Plot

plt.figure(figsize=(5, 4))

data = meta.join(umap)
data['color'] = [common.cmap_tissue[x] for x in data['Tissue']]

data = data.sample(frac=1)


plt.scatter(
    x=data["umap1"],
    y=data["umap2"],
    c=data["color"],
    s=1,
    alpha=0.6,
    edgecolors=None,
    rasterized=True
)

handles = []
labels = []
for tis in common.order_tissue:
    handles.append(Patch(color=common.cmap_tissue[tis]))
    labels.append(common.name_map_tissue[tis])

plt.xticks([])
plt.yticks([])

plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')

for sp in plt.gca().spines.values():
    sp.set_visible(False)

plt.subplots_adjust(right=0.8)

plt.legend(handles, labels,
           loc='upper left',
           bbox_to_anchor=[.92, 1])

#plt.show()
plt.savefig('tissue_umap.svg', dpi=600)

# %% Plot again, but fully vectorized

plt.figure(figsize=(5, 4))

data = meta.join(umap)
data['color'] = [common.cmap_tissue[x] for x in data['Tissue']]

data = data.sample(frac=1)


plt.scatter(
    x=data["umap1"],
    y=data["umap2"],
    c=data["color"],
    s=1,
    alpha=0.6,
    edgecolors=None,
    rasterized=False
)

handles = []
labels = []
for tis in common.order_tissue:
    handles.append(Patch(color=common.cmap_tissue[tis]))
    labels.append(common.name_map_tissue[tis])

plt.xticks([])
plt.yticks([])

plt.xlabel('UMAP-1')
plt.ylabel('UMAP-2')

for sp in plt.gca().spines.values():
    sp.set_visible(False)

plt.subplots_adjust(right=0.8)

plt.legend(handles, labels,
           loc='upper left',
           bbox_to_anchor=[.92, 1])

#plt.show()
plt.savefig('tissue_umap_vector.svg')
