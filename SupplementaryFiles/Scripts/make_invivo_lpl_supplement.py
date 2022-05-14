import pandas as pd
from natsort import natsorted
from excel_utilities import safe_len, series_to_col_width, write_sheet

out_xlsx = snakemake.output['xlsx']

# Cell-meta info

meta = pd.read_table(
    "../10xAll/combined_LPL_zinb/data/meta.txt.gz", index_col=0
)

cluster = pd.read_table(
    "../10xAll/combined_LPL_zinb/cluster_together/clusters.txt", index_col=0
)

umap = pd.read_table(
    "../10xAll/combined_LPL_zinb/umap/umap.txt", index_col=0
).rename({
    'umap1': 'UMAP1',
    'umap2': 'UMAP2'
}, axis=1)


meta_all = pd.concat(
    (meta, cluster, umap),
    axis=1
)


# Fix the cell labels - the -N suffix implies 10x channels
# which we don't want to do
def fix_label(row):
    old_index = row.name
    new_index = old_index.split('-')[0]
    new_index = new_index + '-' + row.Sample
    return new_index


meta_all.index = meta_all.apply(fix_label, axis=1)
meta_all.index.name = 'CellID'

meta_all = meta_all[[
    'Genotype',
    'Batch',
    'Tissue',
    'Sample',
    'Cluster',
    'UMAP1',
    'UMAP2',
]]

# Gather cluster 1 vs. All DE results

cluster_de_all = pd.read_table(
    "../10xAll/combined_LPL_zinb/cluster_together/cluster_de_1vAll_edgeR/cluster_de_results.txt.gz"
)

cluster_de = {}

cluster_ids = natsorted(cluster_de_all.cluster.unique())

cluster_id_map = {x: 'LPL-{}'.format(x[1:]) for x in cluster_ids}

for cid in cluster_ids:

    cluster_de[cid] = (
        cluster_de_all.loc[cluster_de_all.cluster == cid]
        .rename({'gene': 'GeneID'}, axis=1)
        .drop('cluster', axis=1)
        .set_index('GeneID')
        .loc[lambda x: x.FDR < .1]
        .sort_values('FDR')
    )

# Gather within-cluster genotype-de results

cluster_geno_de_all = pd.read_table(
    "../10xAll/combined_LPL_zinb/cluster_together/cluster_de_genotype_edgeR/cluster_de_results.txt.gz"
)

cluster_de_geno = {}

for cid in cluster_ids:

    cluster_de_geno[cid] = (
        cluster_geno_de_all.loc[cluster_geno_de_all.cluster == cid]
        .rename({'gene': 'GeneID'}, axis=1)
        .drop('cluster', axis=1)
        .set_index('GeneID')
        .sort_values('FDR')
    )

    cluster_de_geno[cid]['logFC'] *= -1  #Re-orient so pos values -> higher in Ctrl

# Special sheet for GWAS

df2 = cluster_de_geno['c2'][['GeneSymbol', 'logFC', 'FDR']] \
    .rename({'logFC': 'LPL-2 logFC', 'FDR': 'LPL-2 FDR'}, axis=1) \
    .reset_index()

df9 = cluster_de_geno['c9'][['GeneSymbol', 'logFC', 'FDR']] \
    .rename({'logFC': 'LPL-9 logFC', 'FDR': 'LPL-9 FDR'}, axis=1) \
    .reset_index()

df7 = cluster_de_geno['c7'][['GeneSymbol', 'logFC', 'FDR']] \
    .rename({'logFC': 'LPL-7 logFC', 'FDR': 'LPL-7 FDR'}, axis=1) \
    .reset_index()

gwas_genes = open("../GWAS/genes_delange_implicated.txt").readlines()
gwas_genes = {x.strip().upper() for x in gwas_genes}

gwas_df = df2.merge(df9).merge(df7)

gwas_df['GWAS?'] = [True if x.upper() in gwas_genes else False for x in gwas_df.GeneSymbol]

gwas_df['sortkey'] = gwas_df[['LPL-7 FDR', 'LPL-2 FDR', 'LPL-9 FDR']].min(axis=1)

gwas_df = gwas_df.loc[lambda x: x['sortkey'] < .1]
gwas_df = gwas_df.sort_values(['GWAS?', 'sortkey'], ascending=[False, True])
gwas_df = gwas_df.drop('sortkey', axis=1)
gwas_df = gwas_df.set_index('GeneID')

for cid in cluster_de_geno:
    cluster_de_geno[cid] = (
        cluster_de_geno[cid]
        .loc[lambda x: x.FDR < .1]
    )

# Description page:

Description = """
Column descriptions for all tables

Cells Sheet:

            CellID: Unique identifier for each cell. Consists of 10x barcode and sample ID
            Genotype: 'ctrl' denotes IL23R wt/eGFP, 'KO' denotes IL23R eGFP/eGFP
            Batch: Experimental batch (1 or 2)
            Tissue: Tissue from which the cell was extracted
            Sample: Full unique ID for the cell's sample
            Cluster: LPL-specific cluster assignment
            UMAP1: UMAP axis 1 coordinate
            UMAP2: UMAP axis 2 coordinate

Other sheets:

Each of the other sheets contains the output from an edgeR likelihood ratio test

LPL-X vs. All: Differential expression between cells of cluster LPL-X and the remaining cells.
            Positive logFC values are associated with genes with higher expression in LPL-X
            than in the remaining clusters.

LPL-X Ctrl vs. KO: Differential expression between the Ctrl cells of cluster LPL-X and the KO cells of cluster LPL-X.
            Positive logFC values are associated with genes with higher expression Ctrl cells.

GWAS:       LogFC and FDR are copied from the respective LPL-X Ctrl vs. KO sheets.
            GWAS? represents whether or not the gene is associated with an IBD GWAS loci

Common columns:

            GeneID: Unique Ensembl identifier for the gene being tested
            GeneSymbol: Common symbol used for this gene
            logFC: log2 coefficient of the linear model
            logCPM: log2 counts-per-million average expression
            LR: likelihood ratio from which p-values are derived
            PValue: Associated p-value for this gene's coefficient
            FDR: Benjamini-Hochberg False Discover Rate

"""

description_df = pd.DataFrame({
    'xx': Description.split("\n")[1:]
})


with pd.ExcelWriter(out_xlsx) as writer:
    write_sheet(meta_all, 'Cells', writer)
    for cid in cluster_ids:
        write_sheet(cluster_de[cid], '{} vs. All'.format(cluster_id_map[cid]), writer)

    for cid in cluster_ids:
        write_sheet(
            cluster_de_geno[cid], '{} Ctrl vs. KO'.format(cluster_id_map[cid]), writer)

    write_sheet(gwas_df, 'GWAS', writer)

    description_df.to_excel(
        writer, sheet_name='Description', header=False, index=False)
