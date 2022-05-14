import pandas as pd
from excel_utilities import safe_len, series_to_col_width, write_sheet

out_xlsx = snakemake.output['xlsx']


meta_prefilter = pd.read_table(
    "../10xAll/combined_all_scvi/data/meta.txt.gz", index_col=0
)

cluster_prefilter = pd.read_table(
    "../10xAll/combined_all_scvi/cluster_together/clusters.txt", index_col=0
).rename({
    'Cluster': 'Initial Cluster'
}, axis=1)

umap = pd.read_table(
    "../10xAll/combined_all_main/umap/umap.txt", index_col=0
).rename({
    'umap1': 'UMAP1',
    'umap2': 'UMAP2'
}, axis=1)


cluster = pd.read_table(
    "../10xAll/combined_all_main/cluster_together/clusters.txt", index_col=0
).rename({
    'Cluster': 'Final Cluster'
}, axis=1)

filtered = pd.Series(
    ['no' if x in cluster.index else 'yes' for x in meta_prefilter.index],
    index=meta_prefilter.index,
    name='Filtered?',
)

invivo_meta = pd.concat(
    (meta_prefilter, filtered, cluster_prefilter, cluster, umap),
    axis=1
)

invivo_meta.Sample = [x.replace('ko', 'KO') for x in invivo_meta.Sample]

# Fix the cell labels - the -N suffix implies 10x channels which we don't want to do

def fix_label(row):
    old_index = row.name
    new_index = old_index.split('-')[0]
    new_index = new_index + '-' + row.Sample
    return new_index

invivo_meta.index = invivo_meta.apply(fix_label, axis=1)
invivo_meta.index.name = 'CellID'

# Gather Genotype-De Results

tissue_dir = {
    'SI_LPL': 'combined_SI_LPL_main',
    'Spleen': 'combined_Spleen_main',
    'Colon_LPL': 'combined_Colon_LPL_main',
    'Colon_IEL': 'combined_Colon_IEL_main',
}


geno_de = {}
for tiss in tissue_dir:
    geno_de_tiss = (
        pd.read_table(
            "../10xAll/{}/de_all/de_results_edger.txt.gz".format(tissue_dir[tiss]))
        .set_index('gene')
        .loc[lambda x: x.FDR < .1]
        .sort_values('FDR')
        .drop('cluster', axis=1)
    )
    geno_de_tiss['logFC'] = geno_de_tiss['logFC'] * -1 # So that pos values -> higher in Ctrl
    geno_de_tiss.index.name = 'GeneID'
    geno_de[tiss] = geno_de_tiss

# Gather X vs. Spleen results

de_vs_spleen = pd.read_table(
    "../10xAll/combined_all_main/TissueDE/tissue_de_vs_spleen.txt.gz"
)

lpl_vs_spleen = (
    de_vs_spleen.loc[de_vs_spleen.Contrast == 'LPL']
    .rename({'gene': 'GeneID'}, axis=1)
    .drop('Contrast', axis=1)
    .set_index('GeneID')
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
)
iel_vs_spleen = (
    de_vs_spleen.loc[de_vs_spleen.Contrast == 'IEL']
    .rename({'gene': 'GeneID'}, axis=1)
    .drop('Contrast', axis=1)
    .set_index('GeneID')
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
)

# Gather ANOVA results

de_combined = (
    pd.read_table("../10xAll/combined_all_main/TissueDE/tissue_de2.txt.gz", index_col=0)
    .drop('Contrast', axis=1)
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
)

de_combined.index.name = 'GeneID'

# Description page:

Description = """
Column descriptions for all tables

Cells Sheet:

            CellID: Unique identifier for each cell
            Genotype: 'ctrl' denotes IL23R wt/eGFP, 'KO' denotes IL23R eGFP/eGFP
            Batch: Experimental batch (1 or 2)
            Tissue: Tissue from which the cell was extracted
            Sample: Full unique ID for the cell's sample
            Filtered?: Whether or not the cell was filtered from downstream analysis
            Initial Cluster: Initial clustering based on all cells
            Final Cluster: Final clustering based on only cells that were not filtered
            UMAP1: UMAP axis 1 coordinate
            UMAP2: UMAP axis 2 coordinate


Other sheets:

Each of the other sheets contains the output from an edgeR likelihood ratio test

Spleen Ctrl vs KO: Genotype comparison between cells extracted from the Spleen
            Positive logFC associated with higher expression in Ctrl cells (compared with KO Cells)
Colon LPL Ctrl vs KO: Genotype comparison between cells extracted from the Colon Lamina Propria
            Positive logFC associated with higher expression in Ctrl cells (compared with KO Cells)
Small Intestine LPL Ctrl vs KO: Genotype comparison between cells extracted from the Small Intestine Lamina Propria
            Positive logFC associated with higher expression in Ctrl cells (compared with KO Cells)
Colon IEL Ctrl vs KO: Genotype comparison between cells extracted from the Colon Epithelium
            Positive logFC associated with higher expression in Ctrl cells (compared with KO Cells)
LPL vs Spleen: Comparison between cells taken from the Lamina Propria (Colon and Small Intestine) vs. Spleen
            Positive logFC associated with higher expression in Lamina Propria
IEL vs Spleen: Comparison between cells taken from the Colon Epithelium vs. Spleen
            Positive logFC associated with higher expression in Colon Epithelium
Tissue Combined:  ANOVA-like test for tissue-specific expression
            Tests the inclusion of all tissue coefficients simultaneously
            Reference level - Spleen

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
    write_sheet(invivo_meta, 'Cells', writer)
    write_sheet(geno_de['Spleen'], 'Spleen Ctrl vs KO', writer)
    write_sheet(geno_de['Colon_LPL'], 'Colon LPL Ctrl vs KO', writer)
    write_sheet(geno_de['SI_LPL'], 'Small Intestine LPL Ctrl vs KO', writer)
    write_sheet(geno_de['Colon_IEL'], 'Colon IEL Ctrl vs KO', writer)
    write_sheet(lpl_vs_spleen, 'LPL vs Spleen', writer)
    write_sheet(iel_vs_spleen, 'IEL vs Spleen', writer)
    write_sheet(de_combined, 'Tissue Combined', writer)

    description_df.to_excel(
        writer, sheet_name='Description', header=False, index=False)
