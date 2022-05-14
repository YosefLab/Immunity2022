"""
This script is used to combined the results from two 10x runs together

Data is taken from the individual analysis folders
"""

import pandas as pd

# Combine Data
exp_ctrl_file = "../data_ctrl/expression/expression_sym.txt.gz"
exp_ko_file = "../data_ko/expression/expression_sym.txt.gz"

exp_ctrl = pd.read_table(exp_ctrl_file, index_col=0)
exp_ko = pd.read_table(exp_ko_file, index_col=0)

# Need to replace column labels, label exp_ko columns with -2 instead
exp_ko.columns = [x.replace("-1", "-2") for x in exp_ko.columns]

exp_combined = exp_ctrl.merge(exp_ko, how="outer", left_index=True, right_index=True)
exp_combined[exp_combined.isnull()] = 0.0

# Filter and scale
valid_genes = exp_combined.sum(axis=1) > 10
exp_filtered = exp_combined.loc[valid_genes]

num_umi = exp_combined.sum(axis=0)
exp_scaled = exp_filtered.divide(num_umi, axis=1)*4000

exp_scaled.to_csv("expression_scaled.txt.gz", sep="\t", compression="gzip")

# Write counts combined too
exp_combined.to_csv("expression_sym.txt.gz", sep="\t", compression="gzip")

# Combine QC
qc_ctrl_file = "../data_ctrl/qc/qc.txt.gz"
qc_ko_file = "../data_ko/qc/qc.txt.gz"

qc_ctrl = pd.read_table(qc_ctrl_file, index_col=0)
qc_ko = pd.read_table(qc_ko_file, index_col=0)

qc_ko.index = [x.replace("-1", "-2") for x in qc_ko.index]

qc_combined = pd.concat((qc_ctrl, qc_ko), axis=0)
if "qPC1" in qc_combined.columns:
    qc_combined = qc_combined.drop("qPC1", axis=1)
if "qPC2" in qc_combined.columns:
    qc_combined = qc_combined.drop("qPC2", axis=1)
qc_combined.to_csv("qc.txt.gz", sep="\t", compression="gzip")


# Create meta data

genotype = {x: "ctrl" if "-1" in x else "KO" for x in exp_scaled.columns}
genotype = pd.Series(genotype, name="Genotype")

meta = pd.concat((genotype, ), axis=1)
meta.to_csv("meta.txt.gz", sep="\t", compression="gzip")
