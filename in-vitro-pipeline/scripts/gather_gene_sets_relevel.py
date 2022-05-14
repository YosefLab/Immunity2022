import pandas as pd
import gene_enrich

in_file = snakemake.input['input']
out_file = snakemake.output['out']


res = pd.read_table(in_file)


c1 = res.loc[res.contrast == 'Genotypectrl:GFPeGFP_pos_sub']  # The interaction
cx = res.loc[res.contrast == 'GFP_Ind_sub']  # The indepedent results


c1 = c1.loc[c1.FDR < 0.1]

c1 = c1.set_index('primerid')
cx = cx.set_index('primerid')

# In some cases the logFC cannot be computed because the expression
# is 0 in one group.  In this case, substitute the coefD as a proxy
# indication of effect size

ii = c1.logFC.isnull()
c1.loc[ii, 'logFC'] = c1.loc[ii, 'coefD']

ii = cx.logFC.isnull()
cx.loc[ii, 'logFC'] = cx.loc[ii, 'coefD']

# Create the Geno-Ind set

name = 'Geno-Ind'
genes = cx['logFC']
description = 'Genes differential (GFP+ vs. GFP-) independent of IL23R KO'

gs1 = gene_enrich.GeneSet(
    name=name, genes=genes.index.tolist(),
    description=description,
    values=genes.tolist()
)

# Create the Geno-Dep set

name = 'Geno-Dep'
genes = c1['logFC'] * -1  # For sign consistency before 'relevel'
description = 'Genes differential (GFP+ vs. GFP-) in a IL23R-KO dependent manner'

gs2 = gene_enrich.GeneSet(
    name=name, genes=genes.index.tolist(),
    description=description,
    values=genes.tolist()
)

gene_enrich.write_gmt_file(
    [gs1, gs2], out_file, include_values=True
)
