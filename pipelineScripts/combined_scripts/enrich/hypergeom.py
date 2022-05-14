"""
Performs unisgned gene-set enrichment using the hypergeometric test
"""

from gene_enrich.analyses import gene_set_enrichment
from gene_enrich import load_gene_set_gmt
from statsmodels.sandbox.stats.multicomp import multipletests
from tqdm import tqdm
import pandas as pd


gene_lists_file = snakemake.input["gene_lists"]
enrichment_sets_file = snakemake.input["enrichment_sets"]
background_file = snakemake.input["background"]

out_file = snakemake.output["results"]

gene_lists = load_gene_set_gmt(gene_lists_file)
enrichment_sets = load_gene_set_gmt(enrichment_sets_file)


# Load the background
if background_file.endswith(".feather"):
    import feather
    background_genes = feather.read_dataframe(background_file)["index"].tolist()
else:
    background_genes = pd.read_table(background_file, index_col=0) \
        .index.tolist()

background_genes = [x.upper() for x in background_genes]
background_genes = list(set(background_genes))

print(
    "Loaded {} genes for enrichment background.".format(
        len(background_genes)))

results = []
for gl in tqdm(gene_lists):

    result = gene_set_enrichment(
        gl.genes, all_genes=background_genes,
        gene_sets=enrichment_sets,
        correct_false_neg=False
    )

    reject, pv_corrected, aS, aB = multipletests(
        result['pvalue'].values, alpha=0.05, method='fdr_bh')
    result['FDR'] = pv_corrected
    result['Group'] = gl.name
    result.index.name = 'EnrichmentSet'
    result = result.reset_index()

    results.append(result)

results = pd.concat(results, axis=0)

results.to_csv(out_file, sep="\t")

