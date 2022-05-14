import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tqdm import tqdm

tpm_in = snakemake.input['tpm']
qc_in = snakemake.input['qc']
other_meta_in = snakemake.input['otherMeta']

regress_vars = snakemake.params['regressVars']

tpm_out = snakemake.output['tpm']

out_dir = os.path.dirname(tpm_out)
os.makedirs(out_dir, exist_ok=True)

# qc_in = "data_filtered/qc.txt.gz"

tpm = pd.read_table(tpm_in, index_col=0)
qc = pd.read_table(qc_in, index_col=0)
other_meta = pd.read_table(other_meta_in, index_col=0)

other_meta = other_meta.loc[tpm.columns]
qc = qc.loc[tpm.columns]

meta = other_meta.join(qc).loc[tpm.columns]

tpm_scale = np.log2(tpm + 1).T
tpm_regress = pd.DataFrame(index=tpm_scale.index)

# Regress NGenesDetected out of the rest before doing PCs

# regress_vars = ['NGenesDetected', 'QPC1', 'CC1', 'CC2']
x = meta[regress_vars].values

print("Regressing out covariates from each gene...")
for gene in tqdm(tpm_scale.columns):

    y = tpm_scale[[gene]].values

    model = LinearRegression()
    model.fit(x, y)

    pred_y = model.predict(x) - model.intercept_

    resid_y = y - pred_y

    tpm_regress[gene] = resid_y

tpm_regress = tpm_regress.T
tpm_regress = (2**(tpm_regress))-1
tpm_regress.to_csv(tpm_out, sep="\t", compression="gzip")
