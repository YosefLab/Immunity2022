import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

qc_in = snakemake.input['qc']
qc_out = snakemake.output['qc']

# qc_in = "data_filtered/qc.txt.gz"

qc = pd.read_table(qc_in, index_col=0)

qc_regress = pd.DataFrame(index=qc.index)

# Regress NGenesDetected out of the rest before doing PCs

print("Regression out NGenesDetected from each column...")

for col in qc.columns:
    if col == 'NGenesDetected': continue

    x = qc[['NGenesDetected']].values
    y = qc[[col]].values

    model = LinearRegression()
    model.fit(x, y)

    pred_y = model.predict(x)

    resid_y = pred_y - y

    qc_regress[col] = resid_y

# Now perform PCA
print("Computing PCA...")
qc_scaled = (StandardScaler()
             .fit(qc_regress.values)
             .transform(qc_regress.values))

qc_transformed = (PCA(whiten=True)
                  .fit(qc_scaled)
                  .transform(qc_scaled))

qc['QPC1'] = qc_transformed[:, 0]
qc['QPC2'] = qc_transformed[:, 1]

qc.to_csv(qc_out, sep="\t", compression="gzip")
