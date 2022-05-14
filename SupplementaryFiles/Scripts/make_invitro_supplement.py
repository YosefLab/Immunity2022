import pandas as pd
from excel_utilities import safe_len, series_to_col_width, write_sheet

out_xlsx = snakemake.output['xlsx']

# Gather cell-level data

meta = pd.read_table("../in-vitro-pipeline/Th1_Pipeline/data_filtered/meta.txt.gz", index_col=0)
tsne = pd.read_table("../in-vitro-pipeline/Th1_Pipeline/tsne/tsne.txt", index_col=0)
cc = pd.read_table("../in-vitro-pipeline/Th1_Pipeline/cell_cycle/cc_components.txt.gz", index_col=0)

in_vitro_meta = pd.concat(
    (meta, tsne, cc), axis=1
)

# Rename and re-order columns

in_vitro_meta = in_vitro_meta.rename(
    {
        "Mouse_ID": "MouseID",
        "GFP": "GFP",
        "Stimulus": "Stimulus",
        "Well": "Well",
        "Genotype": "Genotype",
        "tsne1": "Tsne1",
        "tsne2": "Tsne2",
        "CC1": "CellCycle1",
        "CC2": "CellCycle2",
    },
    axis=1,
)

in_vitro_meta.index.name = 'CellID'

# Gather cell-level data: Th17

meta = pd.read_table("../in-vitro-pipeline/Th17_Pipeline/data_filtered/meta.txt.gz", index_col=0)
tsne = pd.read_table("../in-vitro-pipeline/Th17_Pipeline/tsne/tsne.txt", index_col=0)
cc = pd.read_table("../in-vitro-pipeline/Th17_Pipeline/cell_cycle/cc_components.txt.gz", index_col=0)

in_vitro_meta_th17 = pd.concat(
    (meta, tsne, cc), axis=1
)

# Rename and re-order columns

in_vitro_meta_th17 = in_vitro_meta_th17.rename(
    {
        "Mouse_ID": "MouseID",
        "GFP": "GFP",
        "Stimulus": "Stimulus",
        "Well": "Well",
        "Genotype": "Genotype",
        "tsne1": "Tsne1",
        "tsne2": "Tsne2",
        "CC1": "CellCycle1",
        "CC2": "CellCycle2",
    },
    axis=1,
)

in_vitro_meta_th17.index.name = 'CellID'


# Gather GFP+ vs GFP- DE results


data_ctrl = (
    pd.read_table(
        "../in-vitro-pipeline/Th1_Pipeline/DE_Th1_ctrl_noqc/results.txt.gz")
    .query('contrast == "GFPeGFP_pos"')
    .drop('contrast', axis=1)
    .rename({'primerid': 'Gene'}, axis=1)
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
    .set_index('Gene')
)

data_KO = (
    pd.read_table(
        "../in-vitro-pipeline/Th1_Pipeline/DE_Th1_KO_noqc/results.txt.gz")
    .query('contrast == "GFPeGFP_pos"')
    .drop('contrast', axis=1)
    .rename({'primerid': 'Gene'}, axis=1)
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
    .set_index('Gene')
)

# Gather Combined model results

de_res_ix = (
    pd.read_table(
        "../in-vitro-pipeline/Th1_Pipeline/DE_Th1_full2_noqpc_relevel/out.txt.gz")
    .query('contrast == "Genotypectrl:GFPeGFP_pos_sub"')
    .drop('contrast', axis=1)
    .rename({'primerid': 'Gene'}, axis=1)
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
    .set_index('Gene')
)

# %% Gather Th1/Th17 results

data_ctrl_th17 = (
    pd.read_table(
        "../in-vitro-pipeline/Th17_Pipeline/DE_Th17_ctrl_noqc/results.txt.gz")
    .query('contrast == "GFPeGFP_pos"')
    .drop('contrast', axis=1)
    .rename({'primerid': 'Gene'}, axis=1)
    .loc[lambda x: x.FDR < .1]
    .sort_values('FDR')
    .set_index('Gene')
)

# Description page:

Description = """
Column descriptions for all tables

Cells Sheet:

            CellID: Unique identifier for each cell
            MouseID: Biological Sample ID
            GFP: Detection of GFP by FACS
            Stimulus: Cytokine stimulus (12+21+23 = IL-12, IL-21, and IL-23)
            Well: Reaction Well
            Genotype: 'ctrl' denotes IL23R wt/eGFP, 'KO' denotes IL23R eGFP/eGFP
            Tsne1: TSNE axis 1 coordinate
            Tsne2: TSNE axis 2 coordinate
            CellCycle1: Inferred cell-cycle covariate 1
            CellCycle2: Inferred cell-cycle covariate 2

Other sheets:

Each sheet represents the output of MAST in computing a differential expression test
See Methods text for full details

Ctrl GFP+ vs GFP-: GFP coefficient (GFP+ vs GFP-), only Ctrl cells
KO GFP+ vs GFP-: GFP coefficient (GFP+ vs GFP-), only KO cells
Combined Model: Analysis of the interaction term (GFP:Genotype) in a combined model
            with both Ctrl and KO cells. Coefficient tested is for GFP+:GenotypeCtrl

Th1: Cells stimulated in Th1-polarizing conditions (IL12 + IL21 + IL23)
Th17: Cells stimulated in Th17-polarizing conditions (IL1 + IL6 + IL23)

Columns for MAST (common) to all three sheets:

            Gene: Gene for which the test was computed
            coefC: continuous coefficient in the MAST hurdle model
            pvalC: p-value associated with the continuous coefficient
            coefD: discrete coefficient in the MAST hurdle model
            pvalD: p-value associated with the discrete coefficient
            pvalH: overall p-value associated with the MAST hurdle model
            logFC: log2 effect size associated with the tested coefficient
            FDR: Benjamini-Hochberg False Discovery Rate
"""

description_df = pd.DataFrame({
    'xx': Description.split("\n")[1:]
})

# Output the resulting Excel file

with pd.ExcelWriter(out_xlsx) as writer:
    write_sheet(in_vitro_meta, 'Th1 Cells', writer)
    write_sheet(data_ctrl, 'Th1 Ctrl GFP+ vs GFP-', writer)
    write_sheet(data_KO, 'Th1 KO GFP+ vs GFP-', writer)
    write_sheet(de_res_ix, 'Th1 Combined Model', writer)
    write_sheet(in_vitro_meta_th17, 'Th17 Cells', writer)
    write_sheet(data_ctrl_th17, 'Th17 Ctrl GFP+ vs GFP-', writer)
    description_df.to_excel(
        writer, sheet_name='Description', header=False, index=False)
