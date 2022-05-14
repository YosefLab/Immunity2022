import itertools
import gene_enrich
import numpy as np
import pandas as pd
import math
import os
import subprocess

repo_dir = os.path.dirname(subprocess.check_output(['git', 'rev-parse', '--git-dir'])).decode()

# Gather lists and assemble into an Excel for target identification

output = {}
output_file = snakemake.output[0]


# How to score DE tests?

def score_1vall(logFC, FDR):
    if FDR < 0.1 and logFC > 0:
        return logFC
    return 0.0


def score_de(logFC, FDR):
    if FDR < 0.1 and logFC < 0:
        return logFC * -1
    return 0.0


def sign(x):
    return 1 if x >= 0 else -1

# Column 1 - Marker genes cluster c2 of Colon LPL

data = pd.read_table(
    "../combined_LPL/cluster_together/cluster_de_1vAll_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c2'].set_index('GeneSymbol')

# I checked that these are not significant
data = data.loc[~data.index.duplicated()]

output['c2_logFC'] = data.logFC
output['c2_FDR'] = data.FDR
output['c2_Sig'] = pd.Series(
    {gene: score_1vall(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)

# Column 2 - DE Genes within cluster c2 of Colon LPL

data = pd.read_table(
    "../combined_LPL/cluster_together/cluster_de_genotype_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c2'].set_index('GeneSymbol')

# I checked that these are not significant
data = data.loc[~data.index.duplicated()]

output['c2_Geno_logFC'] = data.logFC
output['c2_Geno_FDR'] = data.FDR
output['c2_Geno_Sig'] = pd.Series(
    {gene: score_de(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)
output['c2_Geno_Sign'] = pd.Series(
    {gene: sign(logFC) if sig else 0 for
     gene, logFC, sig in zip(data.index, data.logFC, output['c2_Geno_Sig'])}
)

# Column 2 - Marker genes of cluster c9 of Colon LPL

data = pd.read_table(
    "../combined_LPL/cluster_together/cluster_de_1vAll_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c9'].set_index('GeneSymbol')

# I checked that these are not significant
data = data.loc[~data.index.duplicated()]

output['c9_logFC'] = data.logFC
output['c9_FDR'] = data.FDR
output['c9_Sig'] = pd.Series(
    {gene: score_1vall(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)

# Column 2 - DE Genes within cluster c9 of Colon LPL

data = pd.read_table(
    "../combined_LPL/cluster_together/cluster_de_genotype_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c9'].set_index('GeneSymbol')

# I checked that these are not significant
data = data.loc[~data.index.duplicated()]

output['c9_Geno_logFC'] = data.logFC
output['c9_Geno_FDR'] = data.FDR
output['c9_Geno_Sig'] = pd.Series(
    {gene: score_de(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)
output['c9_Geno_Sign'] = pd.Series(
    {gene: sign(logFC) if sig else 0 for
     gene, logFC, sig in zip(data.index, data.logFC, output['c9_Geno_Sig'])}
)

# Column 2 - Marker genes of cluster c7 of Colon IEL

data = pd.read_table(
    "../combined_Colon_IEL_main/cluster_together/cluster_de_1vAll_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c7'].set_index('GeneSymbol')

# I checked that these are not significant
data = data.loc[~data.index.duplicated()]

output['c7_IEL_logFC'] = data.logFC
output['c7_IEL_FDR'] = data.FDR
output['c7_IEL_Sig'] = pd.Series(
    {gene: score_1vall(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)

# Column 2 - DE Genes within cluster c7 of Colon LPL

data = pd.read_table(
    "../combined_Colon_IEL_main/cluster_together/cluster_de_genotype_edgeR/cluster_de_results.txt.gz")
data = data.loc[data.cluster == 'c7'].set_index('GeneSymbol')

# I checked that these are not significant
data = data.loc[~data.index.duplicated()]

output['c7_IEL_Geno_logFC'] = data.logFC
output['c7_IEL_Geno_FDR'] = data.FDR
output['c7_IEL_Geno_Sig'] = pd.Series(
    {gene: score_de(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)
output['c7_IEL_Geno_Sign'] = pd.Series(
    {gene: sign(logFC) if sig else 0 for
     gene, logFC, sig in
     zip(data.index, data.logFC, output['c7_IEL_Geno_Sig'])}
)

# LPL DE Genes

data = pd.read_table(
    "../combined_LPL/de_all/de_results_edger.txt.gz")
data = data.loc[data.cluster == 'KO_vs_ctrl'].set_index('GeneSymbol')

data = data.loc[~data.index.duplicated()]

output['LPL_Geno_logFC'] = data.logFC
output['LPL_Geno_FDR'] = data.FDR
output['LPL_Geno_Sig'] = pd.Series(
    {gene: score_de(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)
output['LPL_Geno_Sign'] = pd.Series(
    {gene: sign(logFC) if sig else 0 for
     gene, logFC, sig in
     zip(data.index, data.logFC, output['LPL_Geno_Sig'])}
)


# IEL DE Genes

data = pd.read_table(
    "../combined_Colon_IEL_main/de_all/de_results_edger.txt.gz")
data = data.loc[data.cluster == 'KO_vs_ctrl'].set_index('GeneSymbol')

data = data.loc[~data.index.duplicated()]

output['IEL_Geno_logFC'] = data.logFC
output['IEL_Geno_FDR'] = data.FDR
output['IEL_Geno_Sig'] = pd.Series(
    {gene: score_de(logFC, FDR) for
     gene, logFC, FDR in zip(data.index, data.logFC, data.FDR)}
)
output['IEL_Geno_Sign'] = pd.Series(
    {gene: sign(logFC) if sig else 0 for
     gene, logFC, sig in
     zip(data.index, data.logFC, output['IEL_Geno_Sig'])}
)

# Column 3 - In-vitro DE Genes - Th1

gs = gene_enrich.load_gene_set_gmt("../../in-vitro-pipeline/Th1_Pipeline/gene_sets_relevel.gmt")
gs = {x.name: x for x in gs}

output['In-vitro Th1'] = gs['Geno-Ind'].values

# Column 4 - In-vitro DE Genes - Th1 ctrl only

gs = gene_enrich.load_gene_set_gmt("../../in-vitro-pipeline/Th1_Pipeline/gene_sets_relevel.gmt")
gs = {x.name: x for x in gs}

output['In-vitro Th1 IL23R-dependent'] = gs['Geno-Dep'].values * -1

# Column 5/6 - GWAS Genes for IDB (and related)
# Column 5 - 1/0 indicator
# Column 6 - Names of GWAS study(ies)

all_genes = open("../../GWAS/genes_delange_implicated.txt").readlines()
all_genes = [x.strip().upper() for x in all_genes]

gwas_ind = pd.Series({x: 1 for x in all_genes})

output['GWAS?'] = pd.Series(gwas_ind)

# Column 7 TF indicator 0/1

# gs = gene_enrich.load_gene_set_gmt(
#     os.path.join(repo_dir, "Signatures/GO/GO_biological_process_mus_musculus.gmt")
# )
# target_gs = 'regulation of transcription, DNA-templated'

gs = gene_enrich.load_gene_set_gmt(
    os.path.join(repo_dir, "Signatures/GO/GO_molecular_function_mus_musculus.gmt")
)
target_gs = 'transcription factor activity, sequence-specific DNA binding'


gs = [x for x in gs if x.name == target_gs][0]

tf_ind = {x: 1 for x in gs.genes}

output['TF?'] = pd.Series(tf_ind)

# Join it all together in a single table


# First, ensure capitalization isn't an issue
def _fix_series(s):
    s.index = [x.capitalize() for x in s.index]
    s = s[~s.index.duplicated()]
    return s


output = {k: _fix_series(v) for k, v in output.items()}

output_df = pd.DataFrame(output)

output_df.to_csv(output_file, sep='\t')
