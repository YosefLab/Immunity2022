from collections import Counter
import gene_enrich

keywords = {
    'T': ['TH2', 'TFH', 'TH1',
          'TREG', 'TCONV', 'T', 'TH17',
          'TCELL', 'NKTCELL', 'CD4', 'CD8',
          'THYMOCYTE',
          ],

    'Other': ['DC', 'PDC', 'MACROPHAGE',
              'BMDM', 'MONOCYTE', 'NEUTROPHIL', 'MAST', ],

    'B': ['B', 'BCELL', ],
    'NK': ['NKCELL', 'NK', ],
}

all_gs = gene_enrich.load_gene_set_gmt("../C7_IMMSIG_ALL.gmt")


def count_matches(name, terms):
    count = 0
    name_split = name.split("_")
    for term in terms:
        if term in name_split:
            count += 1

    return count


from tqdm import tqdm

gs_counts = {}
for gs in tqdm(all_gs):
    counts = {cat: count_matches(gs.name, terms)
              for cat, terms in keywords.items()}

    gs_counts[gs.name] = counts

# Categorize the gene sets
t_sets = []
b_sets = []
nk_sets = []
other_sets = []
misc_sets = []
cross_sets = []
for gs in gs_counts:
    b = gs_counts[gs]['B']
    t = gs_counts[gs]['T']
    n = gs_counts[gs]['NK']
    other = gs_counts[gs]['Other']

    if b > 0 and t == 0 and n == 0 and other == 0:
        b_sets.append(gs)

    elif b == 0 and t > 0 and n == 0 and other == 0:
        t_sets.append(gs)

    elif b == 0 and t == 0 and n > 0 and other == 0:
        nk_sets.append(gs)

    elif b == 0 and t == 0 and n == 0 and other == 1:  # if other == 2, likely a cross_set
        other_sets.append(gs)

    elif b == 0 and t == 0 and n == 0 and other == 0:
        misc_sets.append(gs)

    else:
        cross_sets.append(gs)

gs_dict = {gs.name: gs for gs in all_gs}


def save_gs(set_names, file_name):
    gene_sets = [gs_dict[gs] for gs in set_names]
    gene_enrich.write_gmt_file(gene_sets, file_name, include_values=True)


save_gs(t_sets, 'TCells.gmt')
save_gs(b_sets, 'BCells.gmt')
save_gs(nk_sets, 'NKCells.gmt')
save_gs(other_sets, 'OtherCells.gmt')
save_gs(misc_sets, 'UnrecognizedCells.gmt')
save_gs(cross_sets, 'CrossSets.gmt')
