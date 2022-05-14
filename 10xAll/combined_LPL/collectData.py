import os
import pandas as pd
from tqdm import tqdm
import shutil
import feather


scaled_file_out = snakemake.output['scaled']
ens_file_out = snakemake.output['ens']
meta_file_out = snakemake.output['meta']
qc_file_out = snakemake.output['qc']
genes_file_out = snakemake.output['genes']

out_dir = os.path.dirname(ens_file_out)
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

ens_to_symbol_file = snakemake.input['genes']
valid_barcodes_file = snakemake.input['valid_barcodes']

valid_barcodes = pd.read_table(valid_barcodes_file, header=None)[0].tolist()
valid_barcodes = pd.Index(valid_barcodes)

data_sets = [
   ('../individual/Colon_IEL_ctrl_2/expression/expression_ens.txt.gz',
    '../individual/Colon_IEL_ctrl_2/qc/qc.txt.gz',
    {'Genotype': 'ctrl', 'Tissue': 'Colon_IEL', 'Batch': '2'}),
   ('../individual/Colon_LPL_ctrl_2/expression/expression_ens.txt.gz',
    '../individual/Colon_LPL_ctrl_2/qc/qc.txt.gz',
    {'Genotype': 'ctrl', 'Tissue': 'Colon_LPL', 'Batch': '2'}),
   ('../individual/SI_LPL_ctrl_2/expression/expression_ens.txt.gz',
    '../individual/SI_LPL_ctrl_2/qc/qc.txt.gz',
    {'Genotype': 'ctrl', 'Tissue': 'SI_LPL', 'Batch': '2'}),
   ('../individual/Spleen_ctrl_2/expression/expression_ens.txt.gz',
    '../individual/Spleen_ctrl_2/qc/qc.txt.gz',
    {'Genotype': 'ctrl', 'Tissue': 'Spleen', 'Batch': '2'}),
   ('../individual/Colon_IEL_KO_2/expression/expression_ens.txt.gz',
    '../individual/Colon_IEL_KO_2/qc/qc.txt.gz',
    {'Genotype': 'KO', 'Tissue': 'Colon_IEL', 'Batch': '2'}),
   ('../individual/Colon_LPL_KO_2/expression/expression_ens.txt.gz',
    '../individual/Colon_LPL_KO_2/qc/qc.txt.gz',
    {'Genotype': 'KO', 'Tissue': 'Colon_LPL', 'Batch': '2'}),
   ('../individual/SI_LPL_KO_2/expression/expression_ens.txt.gz',
    '../individual/SI_LPL_KO_2/qc/qc.txt.gz',
    {'Genotype': 'KO', 'Tissue': 'SI_LPL', 'Batch': '2'}),
   ('../individual/Spleen_KO_2/expression/expression_ens.txt.gz',
    '../individual/Spleen_KO_2/qc/qc.txt.gz',
    {'Genotype': 'KO', 'Tissue': 'Spleen', 'Batch': '2'}),
   ('../individual/Spleen_ctrl_1/expression/expression_ens.txt.gz',
    '../individual/Spleen_ctrl_1/qc/qc.txt.gz',
    {'Genotype': 'ctrl', 'Tissue': 'Spleen', 'Batch': '1'}),
   ('../individual/Spleen_KO_1/expression/expression_ens.txt.gz',
    '../individual/Spleen_KO_1/qc/qc.txt.gz',
    {'Genotype': 'KO', 'Tissue': 'Spleen', 'Batch': '1'}),
   ('../individual/Colon_LPL_ctrl_1/expression/expression_ens.txt.gz',
    '../individual/Colon_LPL_ctrl_1/qc/qc.txt.gz',
    {'Genotype': 'ctrl', 'Tissue': 'Colon_LPL', 'Batch': '1'}),
   ('../individual/Colon_LPL_KO_1/expression/expression_ens.txt.gz',
    '../individual/Colon_LPL_KO_1/qc/qc.txt.gz',
    {'Genotype': 'KO', 'Tissue': 'Colon_LPL', 'Batch': '1'}),
]


all_exp = []
all_meta = []
all_qc = []

for i, ds in tqdm(enumerate(data_sets)):
    index = i+1
    exp_file = ds[0]
    qc_file = ds[1]
    meta_info = ds[2]

    if meta_info['Tissue'] != 'Colon_LPL' and meta_info['Tissue'] != 'SI_LPL':
        continue

    expression = pd.read_table(exp_file, index_col=0)
    qc = pd.read_table(qc_file, index_col=0)

    # replace labels to make each sample unique
    expression.columns = [x.replace("-1", "-"+str(index))
                          for x in expression.columns]
    qc.index = [x.replace("-1", "-"+str(index))
                for x in qc.index]

    # drop qpc if it exists
    if "qPC1" in qc.columns:
        qc = qc.drop("qPC1", axis=1)
    if "qPC2" in qc.columns:
        qc = qc.drop("qPC2", axis=1)

    meta = pd.DataFrame(index=qc.index)

    for var in meta_info:
        meta[var] = meta_info[var]

    # Remove cells not in the 'valid_barcodes' file
    to_keep = expression.columns & valid_barcodes

    expression = expression.loc[:, to_keep]
    qc = qc.loc[to_keep]
    meta = meta.loc[to_keep]

    all_exp.append(expression)
    all_qc.append(qc)
    all_meta.append(meta)


all_exp = pd.concat(all_exp, axis=1)
all_qc = pd.concat(all_qc, axis=0)
all_meta = pd.concat(all_meta, axis=0)

all_meta['Sample'] = [
    '_'.join([a, b, c]) for a, b, c in
    zip(all_meta['Tissue'], all_meta['Genotype'], all_meta['Batch'])
]

# Remove undetected genes
valid_genes = all_exp.sum(axis=1) > 0
all_exp = all_exp.loc[valid_genes]

scaled_exp = all_exp.divide(
    all_exp.sum(axis=0), axis=1) * 3000

# Now, change from ens to symbol
gene_meta = pd.read_table(ens_to_symbol_file, index_col=0, header=None)
gene_meta.columns = ['Symbol']
scaled_exp_sym = scaled_exp.join(gene_meta).groupby('Symbol').sum()

scaled_exp_sym.index.name = 'index'
scaled_exp_sym = scaled_exp_sym.reset_index()
feather.write_dataframe(scaled_exp_sym, scaled_file_out)

all_qc.to_csv(qc_file_out, sep="\t", compression="gzip")
all_meta.to_csv(meta_file_out, sep="\t", compression="gzip")

all_exp.index.name = 'index'
all_exp = all_exp.reset_index()
feather.write_dataframe(all_exp, ens_file_out)

# Copy genes file to the data directory for some other scripts
gene_meta = pd.read_table(ens_to_symbol_file, index_col=0, header=None)
gene_meta.index.name = "geneID"
gene_meta.columns = ["GeneSymbol"]
gene_meta.to_csv(genes_file_out, sep="\t")
