import pandas as pd

from excel_utilities import write_sheet

out_xlsx = snakemake.output['xlsx']

# Load DE results
de = pd.read_table(
    "../cd160_knockout/de_genotype.txt", index_col=0
)
de = de.drop('TDTOMATO', axis=0)  # only existed in one sample
de.index.name = 'GeneSymbol'

Description = """
Column descriptions for all tables

Contains results of the test comparisong CD160 knockout samples vs. Ctrl

Positive LogFC values represent higher expression in the CD160-Knockout group

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

    de = de.loc[lambda x: x['FDR'] < 0.1]
    de = de.sort_values('FDR')
    write_sheet(de, "CD160 KO vs. Ctrl", writer)

    description_df.to_excel(
        writer, sheet_name='Description', header=False, index=False
    )
