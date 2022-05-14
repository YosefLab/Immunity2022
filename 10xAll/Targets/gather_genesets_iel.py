import itertools
import gene_enrich
import pandas as pd
import math

output = snakemake.output['out']

# Gather lists and assemble into an Excel for target identification

out_gs = []

# Set 2 - Marker genes cluster c7 of Colon IEL

data = pd.read_table(
    "../combined_Colon_IEL_main/cluster_together/cluster_de_1vAll_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c7'].set_index('gene')

data = data.loc[data['FDR'] < 0.05]

data = data.sort_values('logFC', ascending=False)

N_GENES = min(200, data.shape[0])

genes = data['GeneSymbol'][:N_GENES].tolist()
genes = list(set(genes))

gs = gene_enrich.GeneSet(
    'c7_IEL', genes, description="Top 200 marker genes for cluster c7 (IEL)"
)

out_gs.append(gs)

# Set 2 - Genotype DE genes of c2 LPL

data = pd.read_table(
    "../combined_Colon_IEL_main/cluster_together/cluster_de_genotype_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c7'].set_index('gene')

data = data.loc[data['FDR'] < 0.05]

data['absFC'] = data['logFC'].abs()

data = data.sort_values('absFC', ascending=False)

N_GENES = min(200, data.shape[0])

genes = data['GeneSymbol'][:N_GENES].tolist()
genes = list(set(genes))

gs = gene_enrich.GeneSet(
    'c7_IEL_geno', genes, description="Top 200 genotype DE genes for cluster c7 (IEL)"
)

out_gs.append(gs)

# Set 3 - Genotype DE genes of LPL (overall)

data = pd.read_table(
    "../combined_Colon_IEL_main/de_all/de_results_edger.txt.gz")
data = data.loc[data.cluster == 'KO_vs_ctrl'].set_index('gene')

data = data.loc[data['FDR'] < 0.05]

data['absFC'] = data['logFC'].abs()

data = data.sort_values('absFC', ascending=False)

N_GENES = min(200, data.shape[0])

genes = data['GeneSymbol'][:N_GENES].tolist()
genes = list(set(genes))

gs = gene_enrich.GeneSet(
    'all_IEL_geno', genes, description="Top 200 genotype DE genes for IEL (all)"
)

out_gs.append(gs)

# Save results

gene_enrich.write_gmt_file(out_gs, output)
