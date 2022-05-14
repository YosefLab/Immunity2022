"""
This script generates how genes were gathered from supplements of de Lange and Huang publications
"""
import pandas as pd
import math

de_lange_supp_table2_file = "../de Lange supplemental tables/NIHMS70681-supplement-Supplementary_Table_2.xlsx"
delange = pd.read_excel(
    de_lange_supp_table2_file,
    skiprows=8, na_values=[], keep_default_na=False,
)


huang_table_file = "../Supplemental Files for Huang et al. Nature 2017/41586_2017_BFnature22969_MOESM2_ESM.xlsx"
huang = pd.read_excel(
    huang_table_file, sheet_name="annotation  summary",
    na_values=[], keep_default_na=False,
)

# Drop a specific locus where there are 89 genes - too ambiguous
ix = delange.loc[
    (delange['Chr'] == 6) &
    (delange['topSNP Position (bp)'] == 32612397)
].index
delange = delange.drop(ix, axis=0)

# Every region in the huang table should be located within a single region in the
# delange table.  Need to assert that this is true and connect the two.

def find_delange_for_huang(row_huang):
    region = row_huang['region (mb)']
    region_chr = int(row_huang['chr'])
    region_start = int(float(region.split('-')[0])*1e6)
    region_end = int(float(region.split('-')[1])*1e6)

    matches = []

    for i, row_delange in delange.iterrows():

        region2_chr = int(row_delange['Chr'])
        region2_loc = int(row_delange['topSNP Position (bp)'])
        region2_in_region = (
            (region2_chr == region_chr) and
            (region2_loc >= region_start) and
            (region2_loc <= region_end)
        )

        if region2_in_region:
            matches.append(i)

    return matches


huang_matches = [
    find_delange_for_huang(row_huang) for _, row_huang in huang.iterrows()
]
huang_matches = [x[0] if len(x) > 0 else -1 for x in huang_matches]

huang['delange_match'] = huang_matches

delange['huang genes'] = ''
for ix, row_huang in huang.iterrows():
    iy = row_huang['delange_match']
    if iy == -1:
        continue
    delange.loc[iy, 'huang genes'] = row_huang['Gene']


# All genes in a delange region
genes_delang_all = set()

# All genes in a delange region except where a gene has been implicated
# (use implicated to override)
genes_delang_implicated = set()

# All genes in a delange region except where a gene has been implicated
# (use implicated to override).  If no gene implicated and region was
# fine-mapped further by Huang, then use the genes in that region instead.
genes_combined = set()

for ix, row in delange.iterrows():
    locus_genes = row['Genes in locus'].split(',')
    locus_genes = [x.strip().upper() for x in locus_genes]
    locus_genes = [x for x in locus_genes if len(x) > 0]

    # Need to get around irregular formatting for this column
    implicated_genes_text = row['Implicated gene']
    implicated_genes = []
    for x in locus_genes:
        if x in implicated_genes_text:
            implicated_genes.append(x)
    if 'ITGB8' in implicated_genes_text:
        implicated_genes.append('ITGB8')  # this one is not listed in the locus genes for some reason


    huang_genes = row['huang genes'].split(',')
    huang_genes = [x.strip().upper() for x in huang_genes]
    huang_genes = [x for x in huang_genes if len(x) > 0]


    genes_delang_all.update(set(locus_genes))

    if len(implicated_genes) > 0:
        genes_delang_implicated.update(set(implicated_genes))
    else:
        genes_delang_implicated.update(set(locus_genes))

    if len(implicated_genes) == 0 and len(huang_genes) > 0:
        genes_combined.update(set(huang_genes))
    elif len(implicated_genes) > 0:
        genes_combined.update(set(implicated_genes))
    else:
        genes_combined.update(set(locus_genes))


with open('genes_delange_all.txt', 'w') as fout:
    for g in genes_delang_all:
        fout.write(g)
        fout.write('\n')

with open('genes_delange_implicated.txt', 'w') as fout:
    for g in genes_delang_implicated:
        fout.write(g)
        fout.write('\n')

with open('genes_delange_w_huang.txt', 'w') as fout:
    for g in genes_combined:
        fout.write(g)
        fout.write('\n')
