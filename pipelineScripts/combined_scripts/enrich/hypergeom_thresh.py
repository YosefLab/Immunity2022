"""
Performs unisgned gene-set enrichment using the hypergeometric test

This version thresholds the background list of genes based on the
lowest-expressed gene in the foreground
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

try:
    out_file_xlsx = snakemake.output['results_xlsx']
except AttributeError:
    out_file_xlsx = None

gene_lists = load_gene_set_gmt(gene_lists_file)
enrichment_sets = load_gene_set_gmt(enrichment_sets_file)


# Load the background
if background_file.endswith(".feather"):
    import feather
    background = feather.read_dataframe(background_file).set_index("index")
else:
    background = pd.read_table(background_file, index_col=0)

bg_mu = background.mean(axis=1)
bg_mu.index = [x.upper() for x in bg_mu.index]


results = []
for gl in tqdm(gene_lists):

    # Compute background
    thresh = bg_mu[list(gl.genes)].min()
    background_genes = bg_mu.index[bg_mu >= thresh]
    background_genes = list(set(background_genes))

    print(
        "{}: Using {} genes for enrichment background.".format(
            gl.name, len(background_genes)))

    # Run enrichment
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

if out_file_xlsx is not None:
    with pd.ExcelWriter(out_file_xlsx) as writer:
        for group in results.Group.unique():
            r_sub = (results.loc[results['Group'] == group]
                     .sort_values('FDR')
                     .drop('Group', axis=1)
                     )
            r_sub = r_sub[
                ['EnrichmentSet', 'fc', 'pvalue', 'FDR', 'matches',
                 'eff_set_size', 'genes']
            ]
            r_sub.to_excel(writer, sheet_name=group, index=False)
