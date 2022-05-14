import numpy as np
import pandas as pd

input_files = snakemake.input
output_file = snakemake.output['out']


# Do TF mapping
# Data downloaded from HGNC - all_genes.txt used for this
tf_gene_map = {
    "Gbx2": "Gbx2", "Gfi1b": "Gfi1b", "Gli1": "Gli1", "Neurod2": "Neurod2",
    "Nr0b1": "Nr0b1", "Prdm14": "Prdm14", "Supt20h": "Supt20h", "Taf7l": "Taf7l",
    "Tp53": "Tp53", "Tp63": "Tp63", "Zfp42": "Zfp42", "Znf143": "Znf143",
    "Znf217": "Znf217", "Znf263": "Znf263", "Znf274": "Znf274", "Znf384": "Znf384",
    "Znf652": "Znf652", "Bp1": "DLX4", "Cbp": "PAG1", "Csb": "ERCC6",
    "E2a": "TCF3", "Eklf": "KLF1", "Er": "EREG", "Eset": "SETDB1",
    "Ews": "EWSR1", "Gf1": "SOS1", "Kap1": "TRIM28", "Nrf2": "NFE2L2",
    "Oct1": "SLC22A1", "P63": "CKAP4", "P68": "POLD3", "Pu.1": "SPI1",
    "Rack7": "ZMYND8", "Ring1b": "RNF2", "Sa1": "STAG1", "Scl": "TAL1",
    "Sfpi1": "SPI1", "Smrt": "NCOR2", "Srebp1": "SREBF1", "Srebp2": "SREBF2",
    "Tbl1": "TBL1Y", "Jarid1a": "KDM5A", "Nfi": "NFIC", "Nmyc": "MYCN",
    "Pcgf4": "BMI1", "Ppar": "PPARA", "Stat5": "STAT5A", "Utx": "KDM6A",
    "Cjun": "JUN", "Cmyc": "MYC", "Egfp-fos": "FOS", "Egfp-gata2": "GATA2",
    "Egfp-hdac8": "Hdac8", "Egfp-junb": "Junb", "Egfp-jund": "Jund", "Egfp-nr4a1": "Nr4a1",
    "Jarid1b-dain": "KDM5B", "Lxr": "NR1H3", "Ncor": "NCOR1", "Nerf2": "ELF2",
    "P300": "EP300", "P53": "TP53", "Pax3-fkhr": "PAX3", "Pbx": "PBX1",
    "Pu": "SPI1", "Pu1": "SPI1", "Smad": "Smad1", "Smad2/3": "Smad2",
    "Tcf12/heb": "Tcf12", "Tcf3/e2a": "Tcf3", "Tcfap2c": "Tfap2c", "Tcfcp2l1": "Tfcp2l1", "Ubf1/2": "Ubtf",
}

tf_gene_map = {x.lower(): y for x, y in tf_gene_map.items()}

def get_tf(esname):
    tf = esname.split('_')[0]
    if tf.lower() in tf_gene_map:
        tf = tf_gene_map[tf.lower()]

    return tf.capitalize()

results = []

for input_file in input_files:
    data = pd.read_table(input_file, index_col=0)
    for group in data.Group.unique():
        data_sub = data.loc[data.Group == group].drop('Group', axis=1)
        data_sub['TF'] = [get_tf(x) for x in data_sub.EnrichmentSet]

        data_sub = (
            data_sub
            .sort_values('FDR')
            .drop_duplicates('TF')
            .set_index('TF')
        )

        # Set the first 20 TF's  with FDR < 0.05 as 'significant'
        data_sub['Sig'] = [x if x > .1 else .1 for x in data_sub['fc']]
        sig1 = data_sub['FDR'] < 0.05
        data_sub.loc[~sig1, 'Sig'] = 0

        # Save columns (TF, FDR, FC, Sig)
        data_sub = data_sub.rename({'fc': 'FC'}, axis=1)

        data_sub = data_sub.loc[:, ['FC', 'FDR', 'Sig']]

        es = 'Encode' if 'encode' in input_file.lower() else 'Chea'
        prefix = '{}_{}'.format(group, es)

        data_sub.columns = ['{}_{}'.format(prefix, x) for x in data_sub.columns]
        results.append(data_sub)

results_df = pd.concat(results, axis=1, sort=True)

genes = set(pd.read_table('targets_pre.txt').iloc[:, 0].tolist())

good_gene = [x for x in results_df.index if x in genes]

results_df = results_df.loc[good_gene]

results_df.to_csv(output_file, sep="\t")
