import pandas as pd
from excel_utilities import safe_len, series_to_col_width, write_sheet

in_file = snakemake.input["table_pre"]
out_xlsx = snakemake.output["xlsx"]


df = pd.read_table(in_file, index_col=0)
df = df.sort_values("Score", ascending=False)

df = df.rename({
    "Score": "Final Score",
    "c9_logFC": "LPL-9_logFC",
    "c9_FDR": "LPL-9_FDR",
    "c9_Geno_logFC": "LPL-9_Geno_logFC",
    "c9_Geno_FDR": "LPL-9_Geno_FDR",
    "c2_logFC": "LPL-2_logFC",
    "c2_FDR": "LPL-2_FDR",
    "c2_Geno_logFC": "LPL-2_Geno_logFC",
    "c2_Geno_FDR": "LPL-2_Geno_FDR",
}, axis=1)

# Re-order columns
df = df[
    [
        "Final Score",
        "In-Vivo Cluster",
        "In-Vivo Tissue",
        "In-Vitro Summary",
        "GWAS?",
        "LPL-9_logFC",
        "LPL-9_FDR",
        "LPL-9_Geno_logFC",
        "LPL-9_Geno_FDR",
        "LPL-2_logFC",
        "LPL-2_FDR",
        "LPL-2_Geno_logFC",
        "LPL-2_Geno_FDR",
        "logFC_LPL_TissueDE",
        "FDR_LPL_TissueDE",
        "LPL_Geno_logFC",
        "LPL_Geno_FDR",
    ]
].copy()

# Description page:

Description = """
Column descriptions for all tables

Targets Sheet:

            Final Score: Total score used to rank genes
            In-Vivo Cluster: Combined score for in-vivo cluster-specific results
            In-Vivo Tissue:  Combined score for in-vivo tissue-specific results
            In-Vitro Summary: Combined score for in-vitro results
            GWAS?: Whether or not the gene is associated with an IBD GWAS loci
            LPL-9_logFC: logFC of the gene's expression in the LPL-9 cluster vs LPL remainder comparison
            LPL-9_FDR: FDR value for the gene in the LPL-9 cluster vs LPL remainder comparison
            LPL-9_Geno_logFC: logFC of the gene's expression in the LPL-9 within-cluster Wild Type vs Knockout comparison
            LPL-9_Geno_FDR: FDR value for the gene in the LPL-9 within-cluster Wild Type vs Knockout comparison
            LPL-2_logFC: logFC of the gene's expression in the LPL-2 cluster vs LPL remainder comparison
            LPL-2_FDR: FDR value for the gene in the LPL-2 cluster vs LPL remainder comparison
            LPL-2_Geno_logFC: logFC of the gene's expression in the LPL-2 within-cluster Wild Type vs Knockout comparison
            LPL-2_Geno_FDR: FDR value for the gene in the LPL-2 within-cluster Wild Type vs Knockout comparison
            logFC_LPL_TissueDE: logFC of the gene's expression in the LPL vs. Spleen comparison
            FDR_LPL_TissueDE: FDR value for the gene in the LPL vs. Spleen comparison
            LPL_Geno_logFC: logFC of the gene's expression in the within-LPL Wild Type vs Knockout comparison
            LPL_Geno_FDR: FDR value for the gene in the within-LPL Wild Type vs Knockout comparison 

All logFC are log2 based.

See methods section for full details.
"""

description_df = pd.DataFrame({
    'xx': Description.split("\n")[1:]
})

with pd.ExcelWriter(out_xlsx) as writer:
    write_sheet(df, 'Targets', writer)

    description_df.to_excel(
        writer, sheet_name='Description', header=False, index=False)
