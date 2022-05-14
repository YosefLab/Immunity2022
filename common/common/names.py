__all__ = ["name_map_tissue", "order_tissue", "order_genotype", "name_map_genotype",
           "name_map_lpl_clusters", "order_lpl_clusters",
           "name_map_iel_clusters", "order_iel_clusters"]
 
name_map_tissue = {
    'Colon_IEL': 'Colon IEL',
    'Colon_LPL': 'Colon LPL',
    'Spleen': 'Spleen',
    'SI_LPL': 'SI LPL',
}

order_tissue = [
    'Spleen',
    'Colon_LPL',
    'SI_LPL',
    'Colon_IEL',
]

name_map_genotype = {
    'ctrl': '$IL23R^{wt/eGFP}$',
    'KO': '$IL23R^{eGFP/eGFP}$',
}

order_genotype = [
    'ctrl', 'KO'
]

order_lpl_clusters = ["c" + str(i) for i in range(12)]

name_map_lpl_clusters = {
    x: 'LPL-{}'.format(x[1:]) for x in order_lpl_clusters
}

order_iel_clusters = ["c" + str(i) for i in range(10)]

name_map_iel_clusters = {
    x: 'IEL-{}'.format(x[1:]) for x in order_iel_clusters
}

