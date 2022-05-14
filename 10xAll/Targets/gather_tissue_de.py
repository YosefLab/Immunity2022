import pandas as pd


in_file = snakemake.input['in_de']
out_file = snakemake.output['out']



#in_file = "../../10x2/combined_all_main/TissueDE/tissue_de.txt.gz"
#out_file ="targets_tissue_pre.txt"

data = pd.read_table(in_file)

data_split = {x: data.loc[data.Contrast == x] for x in data.Contrast.unique()}

def process_set(df, tissue):
    df = df.loc[:, ['GeneSymbol', 'logFC', 'FDR']]
    df = df.loc[~df.GeneSymbol.duplicated()]
    df = df.set_index('GeneSymbol')
    suffix = "_"+tissue+"_TissueDE"
    df.columns = [c+suffix for c in df.columns]
    return df

data_split = {x: process_set(df, x) for x, df in data_split.items()}


data_all = pd.concat(data_split.values(), axis=1, sort=False)


data_all.to_csv(out_file, sep="\t")
