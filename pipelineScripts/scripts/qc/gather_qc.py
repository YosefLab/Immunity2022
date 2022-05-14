# Gathers QC from molecule_stats and Expresion matrix
# Also computes QPCs
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

mol = pd.read_table("../molecule_stats/molecule_qc.txt", index_col=0)
exp = pd.read_table("../expression/expression_sym.txt.gz", index_col=0)


qc = mol


ngenesdetected = (exp > 0).sum(axis=0).rename("num_genes_detected")

qc = pd.concat((qc, ngenesdetected), axis=1)

# Standardize qc
qc_std = qc.subtract(qc.mean(axis=0), axis=1)
qc_std = qc_std.divide(qc_std.std(axis=0), axis=1)

model = PCA(n_components=2)
out = model.fit_transform(qc_std)

out = pd.DataFrame(out, index=qc_std.index, columns=["qPC1", "qPC2"])

qc_all = pd.concat((qc, out), axis=1)

qc_all.to_csv("qc.txt.gz", sep="\t", compression="gzip")
