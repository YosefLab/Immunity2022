import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

tpm_out = snakemake.output['tpm']
counts_out = snakemake.output['counts']
qc_out = snakemake.output['qc']
meta_out = snakemake.output['meta']

# Load data from the two sequencing libraries - FC_02492 and FC_02493

BASE_DIR = '/data/yosef2/Th1_IL23/data'

# FC_02492

COLLECT_DIR = os.path.join(BASE_DIR, 'FC_02492', 'collected', 'rsem')

tpm1 = pd.read_table(os.path.join(COLLECT_DIR, 'rsem_tpmTable_full.txt')) \
    .groupby('Gene Symbol').sum()

counts1 = pd.read_table(
    os.path.join(COLLECT_DIR, 'rsem_readCountsTable_full.txt')) \
    .groupby('Gene Symbol').sum()

qc1 = pd.read_table(os.path.join(COLLECT_DIR, 'qc_full.txt'), index_col=0)

meta1 = pd.read_table(
    os.path.join(COLLECT_DIR, 'sample_meta.txt'),
    index_col=0, dtype=str)

# FC_02493

COLLECT_DIR = os.path.join(BASE_DIR, 'FC_02493', 'collected', 'rsem')

tpm2 = pd.read_table(
    os.path.join(COLLECT_DIR, 'rsem_tpmTable_full.txt')
).groupby('Gene Symbol').sum()

counts2 = pd.read_table(
    os.path.join(COLLECT_DIR, 'rsem_readCountsTable_full.txt')
).groupby('Gene Symbol').sum()

qc2 = pd.read_table(os.path.join(COLLECT_DIR, 'qc_full.txt'), index_col=0)

meta2 = pd.read_table(
    os.path.join(COLLECT_DIR, 'sample_meta.txt'),
    index_col=0, dtype=str)


# Join together
tpm = tpm1.join(tpm2)
counts = counts1.join(counts2)

qc1 = qc1.T
qc2 = qc2.T[qc1.columns]
qc = qc1.append(qc2)

meta2 = meta2[meta1.columns]
meta = meta1.append(meta2)  # append since we are adding rows

# Align all matrices to counts
qc = qc.loc[counts.columns]
meta = meta.loc[counts.columns]
tpm = tpm.loc[:, counts.columns]


# Add a QC Metric
qc['NGenesDetected'] = (tpm > 0).sum(axis=0)

# Remove bad qc metrics
valid_qc = qc.var(axis=0) > 0
qc = qc.loc[:, valid_qc]

# Save results
meta.to_csv(meta_out, sep="\t", compression="gzip")
counts.to_csv(counts_out, sep="\t", compression="gzip")
tpm.to_csv(tpm_out, sep="\t", compression="gzip")
qc.to_csv(qc_out, sep="\t", compression="gzip")
