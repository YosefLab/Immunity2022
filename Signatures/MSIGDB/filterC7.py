import gene_enrich

allc7 = gene_enrich.load_gene_set_gmt("C7_IMMSIG_ALL.gmt")

match_terms = {'TH1', 'TH17', 'TH2', 'CD4', 'CD8', 'Treg', 'teff', 'TCELL', 'NAIVE', 'MEMORY'}
match_terms = {x.lower() for x in match_terms}


def filter(x):
    name_terms = x.name.lower().split("_")
    for i in name_terms:
        if i in match_terms:
            return True

    return False

tcell_c7 = [x for x in allc7 if filter(x)]

gene_enrich.write_gmt_file(tcell_c7, "C7_IMMSIG_TCELL.gmt")
