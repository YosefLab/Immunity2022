# Use this script to extract the signatures from a WikiPathways export
import os
import re
from tqdm import tqdm
import gene_enrich

description = "WikiPathways - Mus Musculus"
gpml_dir = "mus_musculus_wikipathways_export"
out_file = "WikiPathways_Mus_Musculus.gmt"

files_to_scan = os.listdir(gpml_dir)
files_to_scan = [x for x in files_to_scan if x.endswith("gpml")]
files_to_scan = [os.path.join(gpml_dir, x) for x in files_to_scan]

gene_sets = {}

gs_name_re = "Mm_(.+)_WP"
node_re = 'Type="GeneProduct"'
gene_symbol_re = 'TextLabel="([^"]+)"'

for ff in tqdm(files_to_scan):
    mm = re.search(gs_name_re, ff)
    gs_name = mm.group(1)
    gene_sets[gs_name] = []

    for line in open(ff):
        a = re.search(node_re, line)
        if a is None:
            continue

        b = re.search(gene_symbol_re, line)
        if b is None:
            continue

        gene = b.group(1)

        gene_sets[gs_name].append(gene)


gene_set_objs = []
for key in gene_sets:
    genes = list(set(gene_sets[key]))
    if len(genes) == 0:
        continue
    gs = gene_enrich.GeneSet(
        name=key, genes=genes, description=description
    )
    gene_set_objs.append(gs)

gene_enrich.write_gmt_file(
    gene_set_objs, file_name=out_file
)
