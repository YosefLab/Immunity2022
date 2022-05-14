import pandas as pd

from excel_utilities import write_sheet

out_xlsx = snakemake.output['xlsx']

# Load DE results
de_all = pd.read_table(
    "../cytokine_rnaseq/pairwise_de.txt"
)

groups = [
    ("IL-1b+IL-6+IL-23", "TGFb+IL-6"),

    ("TGFb+IL-6", "IL-12"),
    ("IL-1b+IL-6+IL-23", "IL-12"),

    ("TGFb+IL-6", "IL-12+IL-21"),
    ("IL-1b+IL-6+IL-23", "IL-12+IL-21"),

    ("IL-12+IL-21", "IL-12"),
]

Description = """
Column descriptions for all tables

Each of the sheets contains the output from an edgeR likelihood ratio test

The name represents the two groups being compared "e.g., IL-12+IL-21 vs. IL-21"
The order of the name denotes the directionality of the comparison
With the above name, positive LogFC values represent higher expression in the IL-12+IL-21 group

Columns:

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
    for group_num, group_denom in groups:
        de = de_all.loc[
            lambda x: (x['NumComparison'] == group_num) & (x['DenomComparison'] == group_denom)
        ]
        de = de.drop(['NumComparison', 'DenomComparison'], axis=1)
        de = de.loc[lambda x: x['FDR'] < 0.1]
        de = de.sort_values('FDR')
        de = de.set_index('GeneSymbol')
        write_sheet(de, '{} vs {}'.format(group_num, group_denom), writer)

    description_df.to_excel(
        writer, sheet_name='Description', header=False, index=False)
