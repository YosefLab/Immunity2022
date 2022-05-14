import pandas as pd

df1_file = snakemake.input['df1']
df2_file = snakemake.input['df2']
df_de_file = snakemake.input['df_de']

excel_out = snakemake.output['out_excel']
df_out = snakemake.output['out_df']

df1 = pd.read_table(df1_file, index_col=0)
df2 = pd.read_table(df2_file, index_col=0)
df_de = pd.read_table(df_de_file, index_col=0)

for x in (df1, df2, df_de):
    x.index = [y.capitalize() for y in x.index]

# df1 is our gene expression df.  Only retain genes that are in this one
common = df2.index & df1.index
df2 = df2.loc[common]


df_de = df_de.loc[~df_de.index.duplicated()]

# Merge together

df = pd.concat((df1, df2, df_de), axis=1, sort=False)

# Concatenate sig for TF across Encode and Chea
df["c2_LPL_TF_Sig"] = df[[
    "c2_LPL_Chea_Sig",
    "c2_LPL_Encode_Sig",
]].max(axis=1)


df["c7_IEL_TF_Sig"] = df[[
    "c7_IEL_Encode_Sig",
    "c7_IEL_Chea_Sig",
]].max(axis=1)

df["LPL_Geno_TF_Sig"] = df[[
    "all_LPL_geno_Encode_Sig",
    "all_LPL_geno_Chea_Sig",
]].max(axis=1)

df["IEL_Geno_TF_Sig"] = df[[
    "all_IEL_geno_Chea_Sig",
    "all_IEL_geno_Encode_Sig",
]].max(axis=1)


df["In-vitro_Th1_IL23R-dep_TF"] = df[[
    "Geno-Dep_Chea_Sig",
    "Geno-Dep_Encode_Sig",
]].max(axis=1)


df["In-vitro Th1 Sig"] = df["In-vitro Th1"]
df["In-vitro Th1 Dep Sig"] = df["In-vitro Th1 IL23R-dependent"]

def resolve_tissue_de(logfc, fdr):
    if fdr < 0.1 and logfc > 0:
        return logfc
    return 0.0

df["IEL_Tissue_DE"] = [resolve_tissue_de(logfc, fdr) for
                       logfc, fdr in
                       zip(df["logFC_IEL_TissueDE"], df["FDR_IEL_TissueDE"])
                       ]

df["LPL_Tissue_DE"] = [resolve_tissue_de(logfc, fdr) for
                       logfc, fdr in
                       zip(df["logFC_LPL_TissueDE"], df["FDR_LPL_TissueDE"])
                       ]

# df["Spleen_Tissue_DE"] = [resolve_tissue_de(logfc, fdr) for
#                           logfc, fdr in
#                           zip(df["logFC_Spleen_TissueDE"],
#                               df["FDR_Spleen_TissueDE"])
#                           ]


df["Score"] = 0


# Select columns
groups = []

# Sig Summary Columns
groups.append([
    "c2_Sig",
    "c2_Geno_Sig",
    "c2_LPL_TF_Sig",

    "c7_IEL_Sig",
    "c7_IEL_Geno_Sig",
    "c7_IEL_TF_Sig",

    "LPL_Geno_Sig",
    "LPL_Geno_TF_Sig",

    "IEL_Geno_Sig",
    "IEL_Geno_TF_Sig",

    "In-vitro Th1 Sig",
    "In-vitro Th1 Dep Sig",
    "In-vitro_Th1_IL23R-dep_TF",

    "IEL_Tissue_DE",
    "LPL_Tissue_DE",
#     "Spleen_Tissue_DE",

    "GWAS?",
    "TF?",

    "Exp_Sig",

    "Score",

])


# Other info

groups.append([
    "c2_logFC",
    "c2_FDR",
])

groups.append([
    "c2_Geno_logFC",
    "c2_Geno_FDR",
    "c2_Geno_Sign",
])

groups.append([
    "c9_logFC",
    "c9_FDR",
    "c9_Sig",
])

groups.append([
    "c9_Geno_logFC",
    "c9_Geno_FDR",
    "c9_Geno_Sig",
    "c9_Geno_Sign",
])

groups.append([
    "c7_IEL_logFC",
    "c7_IEL_FDR",
])

groups.append([
    "c7_IEL_Geno_logFC",
    "c7_IEL_Geno_FDR",
    "c7_IEL_Geno_Sign",
])

groups.append([
    "LPL_Geno_logFC",
    "LPL_Geno_FDR",
    "LPL_Geno_Sign",
])

groups.append([
    "IEL_Geno_logFC",
    "IEL_Geno_FDR",
    "IEL_Geno_Sign",
])

groups.append([
    "In-vitro Th1",
    "In-vitro Th1 IL23R-dependent",
])

groups.append([
    "c2_LPL_Encode_FC",
    "c2_LPL_Encode_FDR",
    "c2_LPL_Encode_Sig",
    "c2_LPL_Chea_FC",
    "c2_LPL_Chea_FDR",
    "c2_LPL_Chea_Sig",
])

groups.append([
    "c2_LPL_geno_Encode_FC",
    "c2_LPL_geno_Encode_FDR",
    "c2_LPL_geno_Encode_Sig",
    "c2_LPL_geno_Chea_FC",
    "c2_LPL_geno_Chea_FDR",
    "c2_LPL_geno_Chea_Sig",
])

groups.append([
    "c7_IEL_Encode_FC",
    "c7_IEL_Encode_FDR",
    "c7_IEL_Encode_Sig",
    "c7_IEL_Chea_FC",
    "c7_IEL_Chea_FDR",
    "c7_IEL_Chea_Sig",
])

groups.append([
    "c7_IEL_geno_Encode_FC",
    "c7_IEL_geno_Encode_FDR",
    "c7_IEL_geno_Encode_Sig",
    "c7_IEL_geno_Chea_FC",
    "c7_IEL_geno_Chea_FDR",
    "c7_IEL_geno_Chea_Sig",
])

groups.append([
    "all_IEL_geno_Encode_FC",
    "all_IEL_geno_Encode_FDR",
    "all_IEL_geno_Encode_Sig",
    "all_IEL_geno_Chea_FC",
    "all_IEL_geno_Chea_FDR",
    "all_IEL_geno_Chea_Sig",
])

groups.append([
    "all_LPL_geno_Encode_FC",
    "all_LPL_geno_Encode_FDR",
    "all_LPL_geno_Encode_Sig",
    "all_LPL_geno_Chea_FC",
    "all_LPL_geno_Chea_FDR",
    "all_LPL_geno_Chea_Sig",
])


# In-vitro
groups.append([
    "Geno-Dep_Chea_FC",
    "Geno-Dep_Chea_FDR",
    "Geno-Dep_Chea_Sig",
    "Geno-Dep_Encode_FC",
    "Geno-Dep_Encode_FDR",
    "Geno-Dep_Encode_Sig",
])

# Tissue DE
groups.append([
    "logFC_IEL_TissueDE",
    "FDR_IEL_TissueDE",
    "logFC_LPL_TissueDE",
    "FDR_LPL_TissueDE",
#     "logFC_Spleen_TissueDE",
#     "FDR_Spleen_TissueDE",
])

# MegaTable
groups.append([
    "Exp_IL23_4hr_vs_0hr_FDR",
    "Exp_IL23_4hr_vs_0hr_logFC",
    "Exp_IL23_4hr_vs_0hr_Sig",
    "Exp_IL23_20hr_vs_0hr_FDR",
    "Exp_IL23_20hr_vs_0hr_logFC",
    "Exp_IL23_20hr_vs_0hr_Sig",
    "Exp_4hr_IL23+a3/28_vs_a3/28_FDR",
    "Exp_4hr_IL23+a3/28_vs_a3/28_logFC",
    "Exp_4hr_IL23+a3/28_vs_a3/28_Sig",
    "Exp_20hr_IL23+a3/28_vs_a3/28_FDR",
    "Exp_20hr_IL23+a3/28_vs_a3/28_logFC",
    "Exp_20hr_IL23+a3/28_vs_a3/28_Sig",
])

# Zero out all missing values in group 1
# Makes excel calculation easier as min/max will ignore
# The missing values and supply the maximal values erroneously

for col in groups[0]:
    df.loc[df[col].isnull(), col] = 0

# Collect sizes for each group to be used in color formatting later

from itertools import chain
cols = list(chain(*groups))

ranges = []
last_end = 0 
for g in groups:
    r_start = last_end + 1
    r_end = len(g) + r_start - 1
    ranges.append((r_start, r_end))
    last_end = r_end



multi_map = {
    "c2_Sig": ("LPL c2", "1vAll"),
    "c2_Geno_Sig": ("LPL c2", "Geno"),
    "c2_LPL_TF_Sig": ("LPL c2", "TF"),

    "c7_IEL_Sig": ("IEL c7", "1vAll"),
    "c7_IEL_Geno_Sig": ("IEL c7", "Geno"),
    "c7_IEL_TF_Sig": ("IEL c7", "TF"),

    "LPL_Geno_Sig": ("LPL Genotype", "DE"),
    "LPL_Geno_TF_Sig": ("LPL Genotype", "DE_TF"),

    "IEL_Geno_Sig": ("IEL Genotype", "DE"),
    "IEL_Geno_TF_Sig": ("IEL Genotype", "DE_TF"),

    "In-vitro Th1 Sig": ("In-vitro Th1 GFP +/-", "IL23R-indep"),
    "In-vitro Th1 Dep Sig": ("In-vitro Th1 GFP +/-", "IL23R-dep"),
    "In-vitro_Th1_IL23R-dep_TF": ("In-vitro Th1 GFP +/-", "IL23R-dep_TF"),

    "GWAS?": ("Other", "GWAS?"),
    "TF?": ("Other", "TF?"),

    "Score": ("Score", "Score"),

    "c2_logFC": ("LPL c2 1vAll", "logFC"),
    "c2_FDR": ("LPL c2 1vAll", "FDR"),

    "c2_Geno_logFC": ("LPL c2 Geno", "logFC"),
    "c2_Geno_FDR": ("LPL c2 Geno", "FDR"),
    "c2_Geno_Sign": ("LPL c2 Geno", "Sign"),

    "c9_logFC": ("LPL c9 1vAll", "logFC"),
    "c9_FDR": ("LPL c9 1vAll", "FDR"),
    "c9_Sig": ("LPL c9 1vAll", "Sig"),

    "c9_Geno_logFC": ("LPL c9 Geno", "logFC"),
    "c9_Geno_FDR": ("LPL c9 Geno", "FDR"),
    "c9_Geno_Sig": ("LPL c9 Geno", "Sig"),
    "c9_Geno_Sign": ("LPL c9 Geno", "Sign"),

    "c7_IEL_logFC": ("IEL c7 1vAll", "logFC"),
    "c7_IEL_FDR": ("IEL c7 1vAll", "FDR"),
    "c7_IEL_Geno_logFC": ("IEL c7 Geno", "logFC"),
    "c7_IEL_Geno_FDR": ("IEL c7 Geno", "FDR"),
    "c7_IEL_Geno_Sign": ("IEL c7 Geno", "Sign"),

    "LPL_Geno_logFC": ("LPL Geno", "logFC"),
    "LPL_Geno_FDR": ("LPL Geno", "FDR"),
    "LPL_Geno_Sign": ("LPL Geno", "Sign"),
    "IEL_Geno_logFC": ("IEL Geno", "logFC"),
    "IEL_Geno_FDR": ("IEL Geno", "FDR"),
    "IEL_Geno_Sign": ("IEL Geno", "Sign"),

    "c2_LPL_Encode_FC": ("c2_LPL_1vAll_TF", "Encode_FC"),
    "c2_LPL_Encode_FDR": ("c2_LPL_1vAll_TF", "Encode_FDR"),
    "c2_LPL_Encode_Sig": ("c2_LPL_1vAll_TF", "Encode_Sig"),
    "c2_LPL_Chea_FC": ("c2_LPL_1vAll_TF", "Chea_FC"),
    "c2_LPL_Chea_FDR": ("c2_LPL_1vAll_TF", "Chea_FDR"),
    "c2_LPL_Chea_Sig": ("c2_LPL_1vAll_TF", "Chea_Sig"),

    "c2_LPL_geno_Encode_FC": ("c2_LPL_geno_TF", "Encode_FC"),
    "c2_LPL_geno_Encode_FDR": ("c2_LPL_geno_TF", "Encode_FDR"),
    "c2_LPL_geno_Encode_Sig": ("c2_LPL_geno_TF", "Encode_Sig"),
    "c2_LPL_geno_Chea_FC": ("c2_LPL_geno_TF", "Chea_FC"),
    "c2_LPL_geno_Chea_FDR": ("c2_LPL_geno_TF", "Chea_FDR"),
    "c2_LPL_geno_Chea_Sig": ("c2_LPL_geno_TF", "Chea_Sig"),

    "c7_IEL_Encode_FC": ("c7_IEL_1vAll_TF", "Encode_FC"),
    "c7_IEL_Encode_FDR": ("c7_IEL_1vAll_TF", "Encode_FDR"),
    "c7_IEL_Encode_Sig": ("c7_IEL_1vAll_TF", "Encode_Sig"),
    "c7_IEL_Chea_FC": ("c7_IEL_1vAll_TF", "Chea_FC"),
    "c7_IEL_Chea_FDR": ("c7_IEL_1vAll_TF", "Chea_FDR"),
    "c7_IEL_Chea_Sig": ("c7_IEL_1vAll_TF", "Chea_Sig"),

    "c7_IEL_geno_Encode_FC": ("c7_IEL_geno_TF", "Encode_FC"),
    "c7_IEL_geno_Encode_FDR": ("c7_IEL_geno_TF", "Encode_FDR"),
    "c7_IEL_geno_Encode_Sig": ("c7_IEL_geno_TF", "Encode_Sig"),
    "c7_IEL_geno_Chea_FC": ("c7_IEL_geno_TF", "Chea_FC"),
    "c7_IEL_geno_Chea_FDR": ("c7_IEL_geno_TF", "Chea_FDR"),
    "c7_IEL_geno_Chea_Sig": ("c7_IEL_geno_TF", "Chea_Sig"),

    "all_IEL_geno_Encode_FC": ("all_IEL_geno_TF", "Encode_FC"),
    "all_IEL_geno_Encode_FDR": ("all_IEL_geno_TF", "Encode_FDR"),
    "all_IEL_geno_Encode_Sig": ("all_IEL_geno_TF", "Encode_Sig"),
    "all_IEL_geno_Chea_FC": ("all_IEL_geno_TF", "Chea_FC"),
    "all_IEL_geno_Chea_FDR": ("all_IEL_geno_TF", "Chea_FDR"),
    "all_IEL_geno_Chea_Sig": ("all_IEL_geno_TF", "Chea_Sig"),

    "all_LPL_geno_Encode_FC": ("all_LPL_geno_TF", "Encode_FC"),
    "all_LPL_geno_Encode_FDR": ("all_LPL_geno_TF", "Encode_FDR"),
    "all_LPL_geno_Encode_Sig": ("all_LPL_geno_TF", "Encode_Sig"),
    "all_LPL_geno_Chea_FC": ("all_LPL_geno_TF", "Chea_FC"),
    "all_LPL_geno_Chea_FDR": ("all_LPL_geno_TF", "Chea_FDR"),
    "all_LPL_geno_Chea_Sig": ("all_LPL_geno_TF", "Chea_Sig"),


    "In-vitro Th1": ("In-vitro", "IL23R-indep"),
    "In-vitro Th1 IL23R-dependent": ("In-vitro", "IL23R-dep"),

    "Geno-Dep_Chea_FC": ("Geno-Dep_TF", "Chea_FC"),
    "Geno-Dep_Chea_FDR": ("Geno-Dep_TF", "Chea_FDR"),
    "Geno-Dep_Chea_Sig": ("Geno-Dep_TF", "Chea_Sig"),
    "Geno-Dep_Encode_FC": ("Geno-Dep_TF", "Encode_FC"),
    "Geno-Dep_Encode_FDR": ("Geno-Dep_TF", "Encode_FDR"),
    "Geno-Dep_Encode_Sig": ("Geno-Dep_TF", "Encode_Sig"),

    "IEL_Tissue_DE": ("Tissue DE", "IEL"),
    "LPL_Tissue_DE": ("Tissue DE", "LPL"),
#     "Spleen_Tissue_DE": ("Tissue DE", "Spleen"),

    "logFC_LPL_TissueDE": ("Tissue_DE", "LPL logFC"),
    "FDR_LPL_TissueDE": ("Tissue_DE", "LPL FDR"),
    "logFC_IEL_TissueDE": ("Tissue_DE", "IEL logFC"),
    "FDR_IEL_TissueDE": ("Tissue_DE", "IEL FDR"),
#     "logFC_Spleen_TissueDE": ("Tissue_DE", "Spleen logFC"),
#     "FDR_Spleen_TissueDE": ("Tissue_DE", "Spleen FDR"),
}

excel_df = df[cols]

excel_df.columns = (
    pd.MultiIndex.from_tuples(
        [multi_map[x] for x in excel_df.columns]
    )
)


writer = pd.ExcelWriter('Targets.xlsx', engine='xlsxwriter')
excel_df.to_excel(writer, sheet_name='Sheet1')

# Create some formats

light_orange = "#fcd5b4"
light_blue = "#c5d9f1"
light_green = "#d8e4bc"
light_red = "#f2dcdb"

col_format1 = writer.book.add_format({'bg_color': light_blue})
col_format2 = writer.book.add_format({'bg_color': light_green})
col_format3 = writer.book.add_format({'bg_color': light_red})

col_format_n = writer.book.add_format({'bg_color': "#eeeeee"})
col_format_sc = writer.book.add_format({
    'bg_color': '#16365c',
    'font_color': '#ffffff',
    'bold': True,
})

col_formats = [col_format1, col_format2, col_format3]

worksheet = writer.book.worksheets()[0]

for i, r in enumerate(ranges):

    col_format = col_formats[i % len(col_formats)]

    r_start = r[0]
    r_end = r[1]

    worksheet.set_column(r_start, r_end, width=None, cell_format=col_format)

worksheet.set_column(ranges[0][0], ranges[0][1], width=None, cell_format=col_format_n)
worksheet.set_column(ranges[0][1], ranges[0][1], width=None, cell_format=col_format_sc)

writer.save()


# Save regular dataframe too

df.to_csv(df_out, sep="\t", compression="gzip")
