import matplotlib.pyplot as plt
import seaborn as sns

from .names import order_lpl_clusters, order_iel_clusters

__all__ = ["cmap_tissue", "cmap_genotype", "cmap_lpl_clusters", "cmap_iel_clusters"]

# Tissues
colors = sns.cubehelix_palette(4, start=.95, rot=-.7, dark=.3, light=.7)
cmap_tissue = {
    'Spleen': colors[0],
    'Colon_LPL': colors[1],
    'SI_LPL': colors[2],
    'Colon_IEL': colors[3],
}

# Genotype
cmap_genotype = {
    'ctrl': "#8fc6a4",
    'KO': "#c181a8"
}
cmap_genotype['WT'] = cmap_genotype['ctrl']
cmap_genotype['wt'] = cmap_genotype['ctrl']
cmap_genotype['ko'] = cmap_genotype['KO']

# Clusters
lpl_clusters = order_lpl_clusters
palette = sns.husl_palette(12, s=.7, l=.6)
cmap_lpl_clusters = {cl: co for cl, co in zip(lpl_clusters, palette)}

iel_clusters = order_iel_clusters
palette = sns.husl_palette(10, s=.7, l=.6)
cmap_iel_clusters = {cl: co for cl, co in zip(iel_clusters, palette)}
